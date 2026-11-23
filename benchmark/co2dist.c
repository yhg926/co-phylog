/* Co-phylog benchmark tools. GPL-3.0-or-later; see ../LICENSE. */
#include "common.h"
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

static int compare_name(const void *a, const void *b)
{
    return strcmp(*(char * const *)a, *(char * const *)b);
}

static char *full_path(const char *dir, const char *name)
{
    size_t length = strlen(dir) + strlen(name) + 2;
    char *path = allocate(length, 1);
    snprintf(path, length, "%s/%s", dir, name);
    return path;
}

static void check_compatible(const co_file *a, const co_file *b, const char *name)
{
    if (a->half_context != b->half_context || a->min_count != b->min_count || a->min_quality != b->min_quality)
        fail(2, "incompatible half-context/count/quality settings in %s; rebuild with the same settings", name);
}

int main(int argc, char **argv)
{
    if (argc == 2 && (!strcmp(argv[1], "--version") || !strcmp(argv[1], "-V"))) {
        puts("co-phylog co2dist " COPHYLOG_VERSION); return 0;
    }
    if (argc != 2 || argv[1][0] == '-') {
        bool help = argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h"));
        fputs("Usage: co2dist DIRECTORY\nCompare every nonself pair of benchmark .cob files.\n"
              "Output: sample_a, sample_b, co_distance, different_objects, shared_contexts,\n"
              "contexts_a, contexts_b (tab separated). No shared contexts gives NA.\n", help ? stdout : stderr);
        return help ? 0 : 64;
    }
    DIR *dir = opendir(argv[1]);
    if (!dir) fail(2, "cannot open %s: %s", argv[1], strerror(errno));
    size_t count = 0, capacity = 16;
    char **names = allocate(capacity, sizeof(char *));
    struct dirent *item;
    errno = 0;
    while ((item = readdir(dir))) {
        size_t len = strlen(item->d_name);
        if (len < 4 || strcmp(item->d_name + len - 4, ".cob")) continue;
        if (strpbrk(item->d_name, "\t\n\r")) fail(2, "sample filenames cannot contain tabs or newlines");
        if (count == capacity) {
            if (capacity > SIZE_MAX / 2 / sizeof(char *)) fail(2, "too many input files");
            capacity *= 2;
            char **next = realloc(names, capacity * sizeof(char *));
            if (!next) fail(2, "insufficient memory for file list");
            names = next;
        }
        names[count] = strdup(item->d_name);
        if (!names[count]) fail(2, "insufficient memory for file name");
        ++count;
        errno = 0;
    }
    if (errno) fail(2, "cannot finish reading directory");
    closedir(dir);
    if (count < 2) fail(64, "need at least two benchmark .cob files");
    qsort(names, count, sizeof(char *), compare_name);
    /* Validate every input before emitting results, including the final sample. */
    co_file metadata = {0};
    for (size_t i = 0; i < count; ++i) {
        char *path = full_path(argv[1], names[i]);
        co_file co = read_co(path);
        if (i) check_compatible(&metadata, &co, names[i]);
        else metadata = co;
        free(co.tuples); free(path);
    }
    puts("sample_a\tsample_b\tco_distance\tdifferent_objects\tshared_contexts\tcontexts_a\tcontexts_b");
    uint64_t mask = context_mask(metadata.half_context);
    for (size_t a = 0; a + 1 < count; ++a) {
        char *path = full_path(argv[1], names[a]);
        co_file x = read_co(path);
        free(path);
        check_compatible(&metadata, &x, names[a]);
        for (size_t b = a + 1; b < count; ++b) {
            path = full_path(argv[1], names[b]);
            co_file y = read_co(path);
            free(path);
            check_compatible(&metadata, &y, names[b]);
            uint64_t i = 0, j = 0, shared = 0, different = 0;
            while (i < x.size && j < y.size) {
                uint64_t cx = x.tuples[i] & mask, cy = y.tuples[j] & mask;
                if (cx < cy) ++i;
                else if (cx > cy) ++j;
                else {
                    ++shared;
                    different += x.tuples[i] != y.tuples[j];
                    ++i; ++j;
                }
            }
            printf("%s\t%s\t", names[a], names[b]);
            if (shared) printf("%.12g", (double)different / (double)shared);
            else fputs("NA", stdout);
            printf("\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n", different, shared, x.size, y.size);
            free(y.tuples);
        }
        free(x.tuples);
    }
    for (size_t i = 0; i < count; ++i) free(names[i]);
    free(names);
    if (fflush(stdout) || ferror(stdout)) fail(2, "cannot write distance output");
    return 0;
}
