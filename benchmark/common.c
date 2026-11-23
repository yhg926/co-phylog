/* Co-phylog benchmark tools. GPL-3.0-or-later; see ../LICENSE. */
#include "common.h"
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

void fail(int code, const char *format, ...)
{
    va_list args;
    fputs("Error: ", stderr);
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fputc('\n', stderr);
    exit(code);
}

void *allocate(size_t count, size_t size)
{
    if (size && count > SIZE_MAX / size) fail(2, "allocation size overflow");
    void *p = calloc(count ? count : 1, size);
    if (!p) fail(2, "insufficient memory");
    return p;
}

uint32_t parse_uint(const char *value, uint32_t lower, uint32_t upper, const char *name)
{
    char *end;
    errno = 0;
    unsigned long long parsed = strtoull(value, &end, 10);
    if (!*value || *value < '0' || *value > '9' || *end || errno || parsed < lower || parsed > upper)
        fail(64, "%s must be an integer in [%" PRIu32 ", %" PRIu32 "]", name, lower, upper);
    return (uint32_t)parsed;
}

uint64_t tuple_mask(uint32_t half)
{
    return (UINT64_C(1) << (4 * half + 2)) - 1;
}

uint64_t context_mask(uint32_t half)
{
    return tuple_mask(half) ^ (UINT64_C(3) << (2 * half));
}

uint64_t reverse_complement(uint64_t tuple, uint32_t half)
{
    uint64_t result = 0;
    for (uint32_t i = 0; i < 2 * half + 1; ++i) {
        result = (result << 2) | ((tuple & 3) ^ 3);
        tuple >>= 2;
    }
    return result;
}

static void put_uint(FILE *out, uint64_t value, unsigned bytes)
{
    for (unsigned i = 0; i < bytes; ++i) {
        if (fputc((int)(value & 255), out) == EOF) fail(2, "cannot write output");
        value >>= 8;
    }
}

static uint64_t get_uint(FILE *in, unsigned bytes)
{
    uint64_t result = 0;
    for (unsigned i = 0; i < bytes; ++i) {
        int c = fgetc(in);
        if (c == EOF) fail(2, "truncated context file");
        result |= (uint64_t)c << (8 * i);
    }
    return result;
}

void write_co(FILE *out, const co_file *co)
{
    if (fwrite("COBENCH1", 1, 8, out) != 8) fail(2, "cannot write context header");
    put_uint(out, co->half_context, 4);
    put_uint(out, co->min_count, 4);
    put_uint(out, co->min_quality, 4);
    put_uint(out, 0, 4);
    put_uint(out, co->size, 8);
    for (uint64_t i = 0; i < co->size; ++i) put_uint(out, co->tuples[i], 8);
}

co_file read_co(const char *path)
{
    FILE *in = fopen(path, "rb");
    if (!in) fail(2, "cannot open %s: %s", path, strerror(errno));
    char magic[8];
    if (fread(magic, 1, 8, in) != 8 || memcmp(magic, "COBENCH1", 8))
        fail(2, "%s is not a benchmark .cob file; legacy .co files require the historical tools", path);
    co_file co = {0};
    co.half_context = (uint32_t)get_uint(in, 4);
    co.min_count = (uint32_t)get_uint(in, 4);
    co.min_quality = (uint32_t)get_uint(in, 4);
    uint64_t flags = get_uint(in, 4);
    co.size = get_uint(in, 8);
    if (!co.half_context || co.half_context > 15 || !co.min_count || co.min_quality > 93 || flags)
        fail(2, "invalid metadata in %s", path);
    struct stat st;
    if (fstat(fileno(in), &st) || !S_ISREG(st.st_mode) || st.st_size < 32 ||
        co.size != (uint64_t)(st.st_size - 32) / 8 || (st.st_size - 32) % 8 ||
        co.size > SIZE_MAX / sizeof(uint64_t))
        fail(2, "invalid payload length in %s", path);
    co.tuples = allocate((size_t)co.size, sizeof(uint64_t));
    uint64_t mask = context_mask(co.half_context), previous = 0;
    for (uint64_t i = 0; i < co.size; ++i) {
        uint64_t t = co.tuples[i] = get_uint(in, 8), ctx = t & mask;
        if ((t & ~tuple_mask(co.half_context)) ||
            ctx >= (reverse_complement(t, co.half_context) & mask) || (i && ctx <= previous))
            fail(2, "noncanonical, duplicate or unsorted context in %s", path);
        previous = ctx;
    }
    if (ferror(in) || fclose(in)) fail(2, "cannot read %s", path);
    return co;
}
