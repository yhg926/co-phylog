/* Co-phylog benchmark tools. GPL-3.0-or-later; see ../LICENSE. */
#include "common.h"
#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
    uint64_t key_plus_one;
    uint32_t counts[4];
} entry;

typedef struct {
    entry *entries;
    size_t capacity, used;
} table;

static uint64_t mix(uint64_t x)
{
    x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static size_t slot(const table *t, uint64_t key)
{
    size_t i = (size_t)mix(key) & (t->capacity - 1);
    while (t->entries[i].key_plus_one && t->entries[i].key_plus_one != key)
        i = (i + 1) & (t->capacity - 1);
    return i;
}

static void grow(table *t)
{
    if (t->capacity > SIZE_MAX / 2) fail(2, "context table capacity overflow");
    table next = {allocate(t->capacity * 2, sizeof(entry)), t->capacity * 2, t->used};
    for (size_t i = 0; i < t->capacity; ++i)
        if (t->entries[i].key_plus_one)
            next.entries[slot(&next, t->entries[i].key_plus_one)] = t->entries[i];
    free(t->entries);
    *t = next;
}

static void observe(table *t, uint64_t context, unsigned object, uint32_t minimum)
{
    uint64_t key = context + 1;
    size_t i = slot(t, key);
    if (!t->entries[i].key_plus_one) {
        if (t->used >= t->capacity - t->capacity / 4) {
            grow(t);
            i = slot(t, key);
        }
        t->entries[i].key_plus_one = key;
        ++t->used;
    }
    /* Counts above the decision threshold cannot change any retention decision. */
    if (t->entries[i].counts[object] < minimum) ++t->entries[i].counts[object];
}

static int base_code(unsigned char c)
{
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

static uint64_t sort_mask;
static int compare_tuple(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a & sort_mask, y = *(const uint64_t *)b & sort_mask;
    return (x > y) - (x < y);
}

static void usage(FILE *out)
{
    fputs("Usage: fastq2co [options] -o sample.cob reads.fastq[.gz] ...\n"
          "  -k, --half-context INT  Bases on each side of the object [15; 1..15]\n"
          "  -n, --min-count INT     Minimum occurrences per context-object [3]\n"
          "  -Q, --min-quality INT   Minimum Phred+33 quality of every base [20]\n"
          "  -o, --output FILE       Benchmark context file (not legacy .co)\n"
          "  -h, --help              Show help\n"
          "  -V, --version           Show version\n"
          "Use '-' for one standard-input stream. All input files form one sample.\n"
          "No upper-count filter, high-quality single-occurrence shortcut or fixed\n"
          "read-length buffer is used. Ambiguous bases break windows.\n", out);
}

int main(int argc, char **argv)
{
    co_file co = {.half_context=15, .min_count=3, .min_quality=20};
    const char *output = NULL;
    const struct option options[] = {
        {"half-context", required_argument, NULL, 'k'},
        {"min-count", required_argument, NULL, 'n'},
        {"min-quality", required_argument, NULL, 'Q'},
        {"output", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'}, {NULL, 0, NULL, 0}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "k:n:Q:o:hV", options, NULL)) != -1) {
        switch (opt) {
            case 'k': co.half_context = parse_uint(optarg, 1, 15, "half-context"); break;
            case 'n': co.min_count = parse_uint(optarg, 1, UINT32_MAX, "min-count"); break;
            case 'Q': co.min_quality = parse_uint(optarg, 0, 93, "min-quality"); break;
            case 'o': output = optarg; break;
            case 'h': usage(stdout); return 0;
            case 'V': puts("co-phylog fastq2co " COPHYLOG_VERSION); return 0;
            default: usage(stderr); return 64;
        }
    }
    if (!output || !strcmp(output, "-") || optind == argc) { usage(stderr); return 64; }
    struct stat output_stat;
    bool output_exists = stat(output, &output_stat) == 0;
    unsigned stdin_count = 0;
    for (int f = optind; f < argc; ++f) {
        struct stat input_stat;
        if (!strcmp(argv[f], "-")) {
            if (++stdin_count > 1) fail(64, "standard input can be used only once");
            if (fstat(STDIN_FILENO, &input_stat)) fail(2, "cannot inspect standard input");
        } else if (stat(argv[f], &input_stat)) fail(2, "cannot access %s", argv[f]);
        if (output_exists && input_stat.st_dev == output_stat.st_dev && input_stat.st_ino == output_stat.st_ino)
            fail(64, "output must not overwrite an input file");
    }
    table t = {allocate(1024, sizeof(entry)), 1024, 0};
    uint64_t reads = 0, bases = 0, windows = 0, palindromic = 0;
    uint64_t mask = tuple_mask(co.half_context), cmask = context_mask(co.half_context);
    unsigned k = 2 * co.half_context + 1, shift = 2 * (k - 1);
    for (int f = optind; f < argc; ++f) {
        gzFile in = !strcmp(argv[f], "-") ? gzdopen(dup(STDIN_FILENO), "rb") : gzopen(argv[f], "rb");
        if (!in) fail(2, "cannot open %s", argv[f]);
        gzbuffer(in, 1u << 20);
        kseq_t *seq = kseq_init(in);
        if (!seq) fail(2, "cannot initialize FASTQ reader");
        int status;
        while ((status = ks_getc(seq->f)) >= 0) {
            if (status != '@') fail(2, "expected FASTQ header in %s", argv[f]);
            seq->last_char = status;
            if (kseq_read(seq) < 0) fail(2, "incomplete FASTQ record in %s", argv[f]);
            if (!seq->seq.l || !seq->name.l || seq->seq.l != seq->qual.l)
                fail(2, "expected nonempty FASTQ records in %s", argv[f]);
            ++reads; bases += seq->seq.l;
            uint64_t forward = 0, reverse = 0;
            unsigned valid = 0;
            for (size_t p = 0; p < seq->seq.l; ++p) {
                unsigned char q = (unsigned char)seq->qual.s[p];
                if (q < 33 || q > 126) fail(2, "invalid Phred+33 quality in %s", argv[f]);
                int b = base_code((unsigned char)seq->seq.s[p]);
                if (b < 0 || q - 33u < co.min_quality) { forward = reverse = 0; valid = 0; continue; }
                forward = ((forward << 2) | (unsigned)b) & mask;
                reverse = (reverse >> 2) | ((uint64_t)(b ^ 3) << shift);
                if (valid < k) ++valid;
                if (valid < k) continue;
                uint64_t fc = forward & cmask, rc = reverse & cmask;
                /* A self-reverse-complementary context cannot orient its central base. */
                if (fc == rc) { ++palindromic; continue; }
                uint64_t canonical = fc < rc ? forward : reverse;
                observe(&t, canonical & cmask, (unsigned)((canonical >> (2 * co.half_context)) & 3), co.min_count);
                ++windows;
            }
        }
        int zerror;
        gzerror(in, &zerror);
        if (status != -1 || (zerror != Z_OK && zerror != Z_STREAM_END))
            fail(2, "malformed or truncated FASTQ/gzip input: %s", argv[f]);
        kseq_destroy(seq);
        if (gzclose(in) != Z_OK) fail(2, "cannot finish reading %s", argv[f]);
    }
    if (!reads) fail(2, "no FASTQ records found");
    co.tuples = allocate(t.used, sizeof(uint64_t));
    uint64_t conflicts = 0, unsupported = 0;
    for (size_t i = 0; i < t.capacity; ++i) {
        if (!t.entries[i].key_plus_one) continue;
        unsigned supported = 0, object = 0;
        for (unsigned b = 0; b < 4; ++b)
            if (t.entries[i].counts[b] >= co.min_count) { ++supported; object = b; }
        if (supported == 1) co.tuples[co.size++] = (t.entries[i].key_plus_one - 1) | ((uint64_t)object << (2 * co.half_context));
        else if (supported > 1) ++conflicts;
        else ++unsupported;
    }
    free(t.entries);
    sort_mask = cmask;
    qsort(co.tuples, (size_t)co.size, sizeof(uint64_t), compare_tuple);
    size_t length = strlen(output);
    char *temporary = allocate(length + 16, 1);
    snprintf(temporary, length + 16, "%s.tmp.XXXXXX", output);
    int fd = mkstemp(temporary);
    if (fd < 0) fail(2, "cannot create output beside %s: %s", output, strerror(errno));
    FILE *out = fdopen(fd, "wb");
    if (!out) { close(fd); unlink(temporary); fail(2, "cannot open output stream"); }
    write_co(out, &co);
    if (fflush(out) || fsync(fileno(out)) || fclose(out)) {
        unlink(temporary); fail(2, "cannot finish writing %s", output);
    }
    if (rename(temporary, output)) { unlink(temporary); fail(2, "cannot publish %s", output); }
    free(temporary);
    fprintf(stderr, "version=%s half_context=%u min_count=%u min_quality=%u reads=%" PRIu64
            " bases=%" PRIu64 " accepted_windows=%" PRIu64 " palindromic_windows=%" PRIu64
            " retained_contexts=%" PRIu64 " conflicting_contexts=%" PRIu64 " unsupported_contexts=%" PRIu64 "\n",
            COPHYLOG_VERSION, co.half_context, co.min_count, co.min_quality, reads, bases, windows, palindromic,
            co.size, conflicts, unsupported);
    free(co.tuples);
    return 0;
}
