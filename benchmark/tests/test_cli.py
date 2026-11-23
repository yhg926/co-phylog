"""Independent string-based oracle and CLI regression tests on synthetic inputs."""
import collections
import csv
import gzip
import io
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
READER = ROOT / "bin/fastq2co"
DISTANCE = ROOT / "bin/co2dist"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
HEADER = struct.Struct("<8sIIIIQ")


def reverse(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def encode(sequence):
    value = 0
    for base in sequence:
        value = value * 4 + "ACGT".index(base)
    return value


def oracle(records, half=15, minimum=3, quality=20):
    counts = collections.defaultdict(collections.Counter)
    k = 2 * half + 1
    for sequence, scores in records:
        sequence = sequence.upper()
        for start in range(len(sequence) - k + 1):
            window = sequence[start:start + k]
            if set(window) - set("ACGT") or min(scores[start:start + k]) < quality:
                continue
            other = reverse(window)
            forward_ctx = window[:half] + window[half + 1:]
            reverse_ctx = other[:half] + other[half + 1:]
            if forward_ctx == reverse_ctx:
                continue
            chosen = window if forward_ctx < reverse_ctx else other
            counts[chosen[:half] + chosen[half + 1:]][chosen[half]] += 1
    answer = {}
    for ctx, objects in counts.items():
        supported = [obj for obj, count in objects.items() if count >= minimum]
        if len(supported) == 1:
            answer[ctx] = supported[0]
    return answer


def fastq(records):
    return "".join(f"@r{i}\n{sequence}\n+\n{''.join(chr(q + 33) for q in scores)}\n"
                   for i, (sequence, scores) in enumerate(records)).encode()


def load(path):
    payload = path.read_bytes()
    magic, half, minimum, quality, flags, size = HEADER.unpack_from(payload)
    assert magic == b"COBENCH1" and not flags and len(payload) == 32 + 8 * size
    tuples = struct.unpack_from(f"<{size}Q", payload, 32)
    result = {}
    previous = -1
    for value in tuples:
        mask = ((1 << (4 * half + 2)) - 1) ^ (3 << (2 * half))
        assert (value & mask) > previous
        previous = value & mask
        sequence = "".join("ACGT"[(value >> (2 * i)) & 3] for i in reversed(range(2 * half + 1)))
        ctx = sequence[:half] + sequence[half + 1:]
        assert ctx not in result
        result[ctx] = sequence[half]
    return (half, minimum, quality), result


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)
        self.rng = random.Random(20261123)

    def tearDown(self):
        self.temp.cleanup()

    def run_tool(self, args, input_bytes=None, code=0):
        result = subprocess.run([str(a) for a in args], input=input_bytes, capture_output=True, timeout=60)
        self.assertEqual(result.returncode, code, result.stderr.decode(errors="replace"))
        return result

    def sequence(self, length):
        return "".join(self.rng.choice("ACGT") for _ in range(length))

    def sketch(self, records, name="sample", half=15, minimum=3, quality=20):
        source = self.path / (name + ".fq")
        source.write_bytes(fastq(records))
        dest = self.path / (name + ".cob")
        result = self.run_tool([READER, "-k", half, "-n", minimum, "-Q", quality, "-o", dest, source])
        settings, actual = load(dest)
        self.assertEqual(settings, (half, minimum, quality))
        self.assertEqual(actual, oracle(records, half, minimum, quality))
        return dest, actual, result

    def test_defaults_and_help(self):
        seq = self.sequence(80)
        source = self.path / "reads.fq"
        source.write_bytes(fastq([(seq, [20] * 80)] * 3))
        dest = self.path / "default.cob"
        self.run_tool([READER, "-o", dest, source])
        self.assertEqual(load(dest), ((15, 3, 20), oracle([(seq, [20] * 80)] * 3)))
        for tool in (READER, DISTANCE):
            self.assertIn(b"Usage:", self.run_tool([tool, "--help"]).stdout)
            self.assertIn((ROOT / "VERSION").read_text().strip().encode(), self.run_tool([tool, "--version"]).stdout)
            self.run_tool([tool], code=64)

    def test_quality_boundary_at_every_position(self):
        seq = self.sequence(31)
        for pos in range(31):
            scores = [20] * 31
            scores[pos] = 19
            _, actual, _ = self.sketch([(seq, scores)] * 4, name=f"q{pos}")
            self.assertFalse(actual)
        self.assertTrue(self.sketch([(seq, [20] * 31)] * 3)[1])

    def test_quality_settings(self):
        seq = self.sequence(31)
        for threshold in (0, 1, 10, 20, 30, 49, 50, 93):
            for score in {max(0, threshold - 1), threshold, min(93, threshold + 1)}:
                actual = self.sketch([(seq, [score] * 31)] * 3, quality=threshold)[1]
                self.assertEqual(bool(actual), score >= threshold)

    def test_count_and_no_high_quality_bypass(self):
        seq = self.sequence(31)
        for minimum in (1, 2, 3, 10, 20):
            for count in (minimum - 1, minimum, minimum + 1):
                records = [(seq, [93] * 31)] * count + [("N", [20])]
                actual = self.sketch(records, name=f"n{minimum}_{count}", minimum=minimum)[1]
                self.assertEqual(bool(actual), count >= minimum)

    def test_conflicts_after_occurrence_filter(self):
        seq = self.sequence(31)
        alt = seq[:15] + next(b for b in "ACGT" if b != seq[15]) + seq[16:]
        for n_alt in (0, 1, 2, 3, 4):
            records = [(seq, [30] * 31)] * 3 + [(alt, [30] * 31)] * n_alt
            actual = self.sketch(records, name=f"conflict{n_alt}")[1]
            self.assertEqual(bool(actual), n_alt < 3)

    def test_reverse_complement_and_zero_key(self):
        for seq in (self.sequence(150), "A" * 70):
            records = [(seq, [30] * len(seq))] * 2 + [(reverse(seq), [30] * len(seq))]
            self.assertTrue(self.sketch(records)[1])

    def test_palindromic_context(self):
        left = self.sequence(15)
        seq = left + "A" + reverse(left)
        self.assertFalse(self.sketch([(seq, [30] * 31)] * 3)[1])

    def test_random_oracle_and_rehash(self):
        for half in (1, 3, 9, 15):
            for minimum in (1, 3, 10):
                seq = self.sequence(2100)
                scores = [self.rng.choice([19] + [20] * 40) for _ in seq]
                records = [(seq.lower(), scores)] * minimum
                records += [(reverse(seq), list(reversed(scores)))]
                records += [(seq[:100] + "N" + seq[101:], [30] * len(seq))]
                self.sketch(records, name=f"random{half}_{minimum}", half=half, minimum=minimum)

    def test_input_order_invariance(self):
        seq = self.sequence(160)
        alt = seq[:72] + next(b for b in "ACGT" if b != seq[72]) + seq[73:]
        records = [(seq, [40] * len(seq))] * 4 + [(alt, [40] * len(alt))] * 2
        first = self.sketch(records, "first")[0].read_bytes()
        self.rng.shuffle(records)
        self.assertEqual(first, self.sketch(records, "second")[0].read_bytes())

    def test_gzip_stdin_multifile_and_wrapped_fastq(self):
        seq = self.sequence(81)
        records = [(seq, [30] * len(seq))] * 3
        expected = self.sketch(records, "plain")[0].read_bytes()
        compressed = self.path / "reads.fq.gz"
        compressed.write_bytes(gzip.compress(fastq(records)))
        dest = self.path / "other.cob"
        self.run_tool([READER, "-o", dest, compressed])
        self.assertEqual(expected, dest.read_bytes())
        for payload in (fastq(records), gzip.compress(fastq(records))):
            self.run_tool([READER, "-o", dest, "-"], input_bytes=payload)
            self.assertEqual(expected, dest.read_bytes())
        one, two = self.path / "one.fq", self.path / "two.fq"
        one.write_bytes(fastq(records[:1])); two.write_bytes(fastq(records[1:]))
        self.run_tool([READER, "-o", dest, one, two])
        self.assertEqual(expected, dest.read_bytes())
        wrapped = "".join(f"@r{i}\r\n{seq[:40]}\r\n{seq[40:]}\r\n+\r\n{'?' * 40}\r\n{'?' * 41}\r\n" for i in range(3))
        one.write_bytes(wrapped.encode())
        self.run_tool([READER, "-o", dest, one])
        self.assertEqual(expected, dest.read_bytes())

    def test_long_read(self):
        unit, tail = self.sequence(2048), self.sequence(100)
        long_sequence = unit * 1024 + tail
        source = self.path / "long.fq.gz"
        source.write_bytes(gzip.compress(fastq([(long_sequence, [30] * len(long_sequence))])))
        dest = self.path / "long.cob"
        result = self.run_tool([READER, "-n", 1, "-o", dest, source])
        expected = oracle([(unit * 2 + tail, [30] * (len(unit) * 2 + len(tail)))], minimum=1)
        self.assertEqual(load(dest)[1], expected)
        self.assertIn(f"bases={len(long_sequence)} ".encode(), result.stderr)

    def test_distance_oracle(self):
        seq = self.sequence(2000)
        actuals = {}
        for i in range(4):
            changed = list(seq)
            for p in range(80 + i, 2000, 130 + i):
                changed[p] = "ACGT"[("ACGT".index(changed[p]) + 1) % 4]
            text = "".join(changed)
            dest, actual, _ = self.sketch([(text, [30] * len(text))] * 3, name=f"sample{i}")
            actuals[dest.name] = actual
        output = self.run_tool([DISTANCE, self.path]).stdout.decode()
        rows = list(csv.DictReader(io.StringIO(output), delimiter="\t"))
        self.assertEqual(len(rows), 6)
        for row in rows:
            x, y = actuals[row["sample_a"]], actuals[row["sample_b"]]
            shared = x.keys() & y.keys()
            different = sum(x[c] != y[c] for c in shared)
            self.assertEqual(int(row["shared_contexts"]), len(shared))
            self.assertEqual(int(row["different_objects"]), different)
            self.assertAlmostEqual(float(row["co_distance"]), different / len(shared), places=10)

    def test_no_shared_context_is_unavailable(self):
        self.sketch([("N" * 31, [30] * 31)], "empty")
        self.sketch([("A" * 31, [30] * 31)] * 3, "nonempty")
        row = list(csv.DictReader(io.StringIO(self.run_tool([DISTANCE, self.path]).stdout.decode()), delimiter="\t"))[0]
        self.assertEqual(row["shared_contexts"], "0")
        self.assertEqual(row["co_distance"], "NA")

    def test_invalid_arguments(self):
        for option, values in (("-k", (0, 16, -1)), ("-n", (0, -1, "1.2", "3x", "4294967296")), ("-Q", (-1, 94))):
            for value in values:
                self.run_tool([READER, option, value, "-o", self.path / "bad.cob", "-"], code=64)
        self.run_tool([READER, "--unknown"], code=64)
        self.run_tool([DISTANCE, self.path], code=64)

    def test_bad_inputs_do_not_replace_output(self):
        source, dest = self.path / "bad.fq", self.path / "keep.cob"
        dest.write_bytes(b"keep me")
        bad_inputs = (b"", b">fasta\nACGT\n", b"@r\nACGT\n+\n!!!\n", b"@r\nACGT\n+\n!!!!!\n", b"@r\nACGT\n+\n\x7f!!!\n")
        for payload in bad_inputs:
            source.write_bytes(payload)
            self.run_tool([READER, "-o", dest, source], code=2)
            self.assertEqual(dest.read_bytes(), b"keep me")
        complete = fastq([("A" * 50, [30] * 50)] * 3)
        for suffix in (b"@", b"@incomplete", b"garbage\n", b">fasta\nACTG\n"):
            source.write_bytes(complete + suffix)
            self.run_tool([READER, "-o", dest, source], code=2)
            self.assertEqual(dest.read_bytes(), b"keep me")
        source.write_bytes(gzip.compress(fastq([("A" * 50, [30] * 50)] * 3))[:-5])
        self.run_tool([READER, "-o", dest, source], code=2)
        self.assertEqual(dest.read_bytes(), b"keep me")
        self.run_tool([READER, "-o", dest, self.path / "missing.fq"], code=2)

    def test_inplace_output_rejected(self):
        source = self.path / "reads.fq"
        source.write_bytes(fastq([("A" * 40, [30] * 40)] * 3))
        original = source.read_bytes()
        self.run_tool([READER, "-o", source, source], code=64)
        alias = self.path / "alias.cob"
        alias.symlink_to(source)
        self.run_tool([READER, "-o", alias, source], code=64)
        self.assertEqual(original, source.read_bytes())

    def test_mismatched_metadata_and_corrupt_files(self):
        records = [("A" * 40, [30] * 40)] * 4
        first = self.sketch(records, "first")[0]
        second = self.sketch(records, "second", minimum=4)[0]
        self.assertFalse(self.run_tool([DISTANCE, self.path], code=2).stdout)
        valid = first.read_bytes()
        corruptions = [b"legacy!!", valid[:-1], valid + b"x", valid[:24] + struct.pack("<Q", 2**63) + valid[32:]]
        for corrupt in corruptions:
            second.write_bytes(corrupt)
            self.assertFalse(self.run_tool([DISTANCE, self.path], code=2).stdout)


if __name__ == "__main__":
    unittest.main()
