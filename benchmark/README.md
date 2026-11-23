# Co-phylog benchmark tools

Version 0.1.0-benchmark.1 is an opt-in benchmark implementation of co-phylog's
context-object comparison. It adds explicit quality and occurrence controls
and dynamically sized FASTQ input. The historical tools in the parent
directory are unchanged. This is a prerelease, not a drop-in replacement for
historical binaries or a claim of improved biological accuracy.

Please cite the original method:
Huiguang Yi and Li Jin. *Co-phylog: an assembly-free phylogenomic approach for
closely related organisms.* Nucleic Acids Research 41(7), e75 (2013).
[doi:10.1093/nar/gkt003](https://doi.org/10.1093/nar/gkt003).

## Build and test

Requirements: a C11 compiler, make, zlib development headers, and Python 3 for
tests. On Debian/Ubuntu the corresponding packages are `build-essential`,
`zlib1g-dev`, and `python3`.

From the repository root:

```bash
make -C benchmark
make -C benchmark test
benchmark/bin/fastq2co --version
benchmark/bin/co2dist --version
```

`make -C benchmark sanitize` runs the same tests with AddressSanitizer and
UndefinedBehaviorSanitizer on Linux. LeakSanitizer requires an environment
without ptrace restrictions. Release builds do not require sanitizer runtimes.
External `CC`, `CFLAGS`, `CPPFLAGS`, `LDFLAGS` and `LDLIBS` are respected.

## Quick start

Build one sample from paired FASTQ files, and another from long reads:

```bash
mkdir -p contexts
benchmark/bin/fastq2co --half-context 15 --min-count 3 --min-quality 20 \
  -o contexts/sample_a.cob sample_a_R1.fastq.gz sample_a_R2.fastq.gz
benchmark/bin/fastq2co --half-context 15 --min-count 3 --min-quality 20 \
  -o contexts/sample_b.cob sample_b.fastq.gz
benchmark/bin/co2dist contexts > distances.tsv
```

All input files on one reader command are combined as **one sample**. Read
pairs are not merged or overlap-deduplicated. Plain FASTQ, gzip FASTQ, wrapped
FASTQ, CRLF line endings and standard input are supported:

```bash
gzip -dc sample_a.fastq.gz | benchmark/bin/fastq2co -o contexts/sample_a.cob -
```

To change the minimum occurrence count, use `--min-count 10` or `--min-count 20`
(short form `-n`). Use `--min-quality` (short form `-Q`) to change the quality
threshold. A value of zero disables quality-based rejection of valid Phred+33
bases, but input encoding is still validated.

## Exact retention rules

1. `--half-context 15` means 15 context bases on each side of one central object:
   a **31-base window**. The option accepts 1 through 15; the default is 15.
   Its short form `-k` also specifies **half-context**, not total window size.
2. Every base of a window, including the central object, must have Phred+33
   quality at least `--min-quality` (default **20**, inclusive). A low-quality
   or ambiguous base breaks windows; sequence is never concatenated across it.
3. Orient a window by the lexicographically smaller forward/reverse-complement
   context, ignoring its central object. A context equal to its own reverse
   complement is excluded because it cannot orient the central base uniquely.
4. Count each central base separately within its canonical context. The minimum
   is `--min-count` (default **3**). Counts refer to window observations, not
   distinct reads or independent molecules. Counts saturate at the threshold
   internally; this does not change retention decisions.
5. Keep a context only when **exactly one** central base reaches the threshold.
   An alternative base below the threshold does not cause a conflict. If two
   or more bases reach it, discard the context.

There is no upper-abundance cutoff, inferred genome coverage interval,
high-quality single-occurrence shortcut, sketch subsampling or quality-based
read-length truncation. FASTQ buffers grow dynamically; resource and parser
integer limits still apply. A synthetic 2,097,252-base read is covered by tests.
Memory for counting scales with distinct contexts, not total read yield.

These rules intentionally differ from the historical reader's hard-coded
quality shortcuts and support-bit logic. The benchmark label and version must
be reported when publishing results; do not relabel old outputs as this version.

## Distance output

`co2dist` compares every unordered nonself pair in filename order and writes:

| Column | Definition |
|---|---|
| `sample_a`, `sample_b` | Input filenames, including `.cob` |
| `different_objects` | Shared contexts whose central bases differ |
| `shared_contexts` | Contexts retained by both samples |
| `co_distance` | `different_objects / shared_contexts`, or `NA` when none are shared |
| `contexts_a`, `contexts_b` | Retained contexts in each input |

The numerator is a count of differing central objects on shared contexts,
**not automatically a complete whole-genome SNP count**. Missing contexts and
different quality/callable regions affect which differences can be observed.
No identity calibration or SNP-count calibration is applied.

All files in one comparison must have identical half-context, minimum count
and minimum quality settings. A mismatch is an error. Counting output and
quality parameters are printed to stderr; record the complete command,
version, input checksums and these messages for reproducibility.

## Output safety and format

The `.cob` format is deliberately distinct from historical headerless `.co`
files. Do not pass `.cob` files to the historical programs. The benchmark
distance reader rejects unrecognized formats, invalid lengths, noncanonical
contexts, duplicates, unsorted records and incompatible metadata.

All integer fields are little-endian. The 32-byte header is:

| Bytes | Field |
|---|---|
| 0-7 | ASCII `COBENCH1` |
| 8-11 | Half-context, uint32 |
| 12-15 | Minimum count, uint32 |
| 16-19 | Minimum quality, uint32 |
| 20-23 | Reserved flags, zero |
| 24-31 | Number of records, uint64 |

Each record is an eight-byte 2-bit encoded canonical window (A=0, C=1, G=2,
T=3). Records are sorted by the context, with the central two bits ignored for
ordering. There is exactly one record per retained context, including context
zero when valid. The format carries no base-position coordinate or abundance.

Reader output is written to a temporary file beside the destination and renamed
only after successful parsing and writing. Existing output is replaced on
success. Using an input file as the output, including through a hard link or
symlink, is rejected. Failed writes may leave an identifiable `.tmp.*` file;
these are never used as distance inputs. Empty retained context sets are valid
and produce unavailable distances, not artificial zero distances.

Exit statuses are 0 for success/help/version, 64 for invalid command usage, and
2 for input, format, memory or I/O failures. For piped shell workflows use
`set -euo pipefail`.

## Validation and licensing

Tests compare native results with an independent Python string-based oracle,
including count and quality boundaries, reverse complements, conflicts,
palindromic contexts, context zero, rehashing, input order, long reads,
gzip/stdin/multiple inputs, malformed records, corrupt binary files and no-call
behavior. These tests validate software semantics; they are not a biological
accuracy benchmark.

The benchmark source uses the repository's GPL-3.0-or-later license. The
vendored FASTQ parser retains its MIT license; see [vendor/README.md](vendor/README.md).
