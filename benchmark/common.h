/* Co-phylog benchmark tools. GPL-3.0-or-later; see ../LICENSE. */
#ifndef COPHYLOG_BENCH_COMMON_H
#define COPHYLOG_BENCH_COMMON_H

#include <inttypes.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    uint32_t half_context, min_count, min_quality;
    uint64_t size;
    uint64_t *tuples;
} co_file;

void fail(int code, const char *format, ...);
void *allocate(size_t count, size_t size);
uint32_t parse_uint(const char *value, uint32_t lower, uint32_t upper, const char *name);
uint64_t tuple_mask(uint32_t half);
uint64_t context_mask(uint32_t half);
uint64_t reverse_complement(uint64_t tuple, uint32_t half);
void write_co(FILE *out, const co_file *co);
co_file read_co(const char *path);

#endif
