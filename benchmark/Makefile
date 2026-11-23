CC ?= cc
CFLAGS ?= -O3
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?=
VERSION := $(shell cat VERSION)
WARNINGS = -std=c11 -Wall -Wextra -Wpedantic
DEFINES = -D_POSIX_C_SOURCE=200809L -DCOPHYLOG_VERSION='"$(VERSION)"'

.PHONY: all test clean sanitize
all: bin/fastq2co bin/co2dist

bin/fastq2co: fastq2co.c common.c common.h vendor/kseq.h VERSION | bin
	$(CC) $(CPPFLAGS) $(DEFINES) $(CFLAGS) $(WARNINGS) fastq2co.c common.c $(LDFLAGS) $(LDLIBS) -lz -o $@

bin/co2dist: co2dist.c common.c common.h VERSION | bin
	$(CC) $(CPPFLAGS) $(DEFINES) $(CFLAGS) $(WARNINGS) co2dist.c common.c $(LDFLAGS) $(LDLIBS) -o $@

bin:
	mkdir -p bin

test: all
	python3 -m unittest discover -s tests -v

sanitize:
	$(MAKE) clean
	$(MAKE) CFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie' LDFLAGS='-fsanitize=address,undefined -no-pie' test

clean:
	rm -rf bin
