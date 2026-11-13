"""Generate the bundled demo datasets for the web UI.

Four datasets are produced, each with the four lanes
(R1/R2/I1/I2 .fastq.gz) and a samples.csv:

* normal    - every read matches a unique sample (perfect or 1 mismatch)
* ambiguous - one read with no candidate and one read whose best total
              distance is shared by two samples (a tie); both go to
              undetermined
* corrupt   - same content as normal but the R2 gzip stream is truncated
              at its tail, so the CRC / ISIZE check fails and the whole
              batch fails.
* trimming  - reads whose R1/R2 carry known 3' library adapters, including
              full adapters, an internal adapter, an overtrimmed short pair
              and an undetermined read, so the cutadapt trimming page can be
              exercised end to end.

The sample sheet is chosen so that sampleA and sampleB differ at exactly
two positions on both index ends; flipping one of those positions yields a
read that sits at total distance 2 from each of them -> a genuine tie.

Run:  python -I scripts/generate_samples.py [output_root]
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.join(HERE, "..", "samples")

SAMPLE_SHEET = [
    ("sampleA", "AACGTGAT", "AACGTGAT"),
    ("sampleB", "ATCGTGAA", "ATCGTGAA"),  # 2 flips vs A on each end
    ("sampleC", "GCTACGAA", "TTACGCAG"),
]

# Known 3' library adapters used by the trimming demo.  Submit these on the
# trimming page together with min_length=10.
TRIM_ADAPTER_R1 = "AGATCGGAAGAG"          # 12 bp, 3' adapter of R1
TRIM_ADAPTER_R2 = "CTGTCTCTTATA"          # 12 bp, 3' adapter of R2


def _record(name: str, seq: str) -> bytes:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode("ascii")


def _write_lane(path: str, records: list[bytes]) -> None:
    with gzip.open(path, "wb", compresslevel=6) as fh:
        for rec in records:
            fh.write(rec)


def _read_names(n: int) -> list[str]:
    # Whitespace suffixes exercise the "split at first whitespace" rule.
    return [f"SEQ-{i:04d} run=L1" for i in range(1, n + 1)]


def _build_lanes(
    index_pairs: list[tuple[str, str]], names: list[str], out_dir: str
) -> None:
    r1, r2, i1, i2 = [], [], [], []
    for idx, (name, (i7, i5)) in enumerate(zip(names, index_pairs), start=1):
        base = name.split(" ")[0], name.split(" ", 1)[1]
        # R1/R2 carry /1 and /2 suffixes to confirm normalisation.
        r1.append(_record(f"{base[0]}/1 {base[1]}", f"R1SEQ{idx:04d}AAAA"))
        r2.append(_record(f"{base[0]}/2 {base[1]}", f"R2SEQ{idx:04d}TTTT"))
        i1.append(_record(name, i7))
        i2.append(_record(name, i5))
    _write_lane(os.path.join(out_dir, "R1.fastq.gz"), r1)
    _write_lane(os.path.join(out_dir, "R2.fastq.gz"), r2)
    _write_lane(os.path.join(out_dir, "I1.fastq.gz"), i1)
    _write_lane(os.path.join(out_dir, "I2.fastq.gz"), i2)


def _write_sheet(out_dir: str) -> None:
    with open(os.path.join(out_dir, "samples.csv"), "w", encoding="utf-8") as fh:
        fh.write("sample_id,i7,i5\n")
        for row in SAMPLE_SHEET:
            fh.write(",".join(row) + "\n")


def normal_dataset(out_dir: str) -> list[tuple[str, str]]:
    a7, a5 = SAMPLE_SHEET[0][1], SAMPLE_SHEET[0][2]
    b7, b5 = SAMPLE_SHEET[1][1], SAMPLE_SHEET[1][2]
    c7, c5 = SAMPLE_SHEET[2][1], SAMPLE_SHEET[2][2]
    pairs = [
        (a7, a5),                       # perfect -> sampleA
        ("AACGTGTT", a5),               # 1 mismatch on i7 -> sampleA
        (b7, "ATCGTGAC"),               # 1 mismatch on i5 -> sampleB
        ("ATCGTGAG", "ATCGTGAG"),       # 1+1 mismatch -> sampleB (total 2, unique)
        (c7, c5),                       # perfect -> sampleC
        (a7, a5),                       # perfect -> sampleA
    ]
    _build_lanes(pairs, _read_names(len(pairs)), out_dir)
    _write_sheet(out_dir)
    return pairs


def ambiguous_dataset(out_dir: str) -> None:
    a7, a5 = SAMPLE_SHEET[0][1], SAMPLE_SHEET[0][2]
    b7, b5 = SAMPLE_SHEET[1][1], SAMPLE_SHEET[1][2]
    # A and B differ at positions 1 and 7 on each end.  A read that keeps
    # A's base at position 1 and B's base at position 7 (AACGTGAA) is
    # distance 1 from each on each end -> both pass the <=1-per-end gate
    # and share total distance 2: a genuine tie -> undetermined.
    tie_pair_seq = "AACGTGAA"
    pairs = [
        (a7, a5),                       # unique -> sampleA
        (b7, b5),                       # unique -> sampleB
        (SAMPLE_SHEET[2][1], SAMPLE_SHEET[2][2]),  # unique -> sampleC
        ("GGGGGGGG", "CCCCCCCC"),       # >1 away from every sample -> no candidate
        (tie_pair_seq, tie_pair_seq),   # tied best (2 vs 2) -> undetermined
    ]
    _build_lanes(pairs, _read_names(len(pairs)), out_dir)
    _write_sheet(out_dir)


def corrupt_dataset(out_dir: str, normal_dir: str) -> None:
    # Copy normal, then cut the R2 gzip stream nine bytes short, removing
    # its 8-byte trailer (CRC32 + ISIZE) and one byte of deflate data.
    for name in ("R1.fastq.gz", "I1.fastq.gz", "I2.fastq.gz", "samples.csv"):
        shutil.copyfile(os.path.join(normal_dir, name), os.path.join(out_dir, name))
    src = os.path.join(normal_dir, "R2.fastq.gz")
    dst = os.path.join(out_dir, "R2.fastq.gz")
    with open(src, "rb") as fh:
        data = fh.read()
    with open(dst, "wb") as fh:
        fh.write(data[:-9])


def trimming_dataset(out_dir: str) -> None:
    """Build a batch whose R1/R2 reads embed known 3' adapters.

    Submit adapters ``TRIM_ADAPTER_R1`` / ``TRIM_ADAPTER_R2`` with
    ``min_length=10`` on the trimming page:

    * SEQ-0001 (sampleA): no adapter, full 20 bp kept
    * SEQ-0002 (sampleA): full adapters on both ends, 12 bp insert kept
    * SEQ-0003 (sampleB): 8 bp biological bases then adapters -> trimming
      leaves 8 bp, below the minimum, the whole pair is filtered out
    * SEQ-0004 (sampleB): adapter sequence appears inside the read; the 3'
      suffix from the adapter onward is removed
    * SEQ-0005 (sampleC): R1 has its adapter, R2 has none
    * SEQ-0006 (sampleC): partial adapter ligation (6 of 12 bp) which the
      default overlap rule still recognises
    * SEQ-0007: indices match no sample -> undetermined, but its reads are
      still trimmable
    """
    a7, a5 = SAMPLE_SHEET[0][1], SAMPLE_SHEET[0][2]
    b7, b5 = SAMPLE_SHEET[1][1], SAMPLE_SHEET[1][2]
    c7, c5 = SAMPLE_SHEET[2][1], SAMPLE_SHEET[2][2]
    index_pairs = [
        (a7, a5),                                  # 1 -> sampleA
        (a7, a5),                                  # 2 -> sampleA
        (b7, b5),                                  # 3 -> sampleB (filtered)
        (b7, b5),                                  # 4 -> sampleB
        (c7, c5),                                  # 5 -> sampleC
        (c7, c5),                                  # 6 -> sampleC
        ("GGGGGGGG", "CCCCCCCC"),                  # 7 -> undetermined
    ]
    names = _read_names(len(index_pairs))
    r1, r2, i1, i2 = [], [], [], []

    def rec(name: str, seq: str) -> bytes:
        return _record(name, seq)

    reads = [
        # (R1 sequence, R2 sequence)
        ("AAAACCCCGGGGTTTTAA", "TTTTGGGGCCCCAAAATT"),                       # 1 clean
        ("AAAACCCCAAAA" + TRIM_ADAPTER_R1,
         "TTTTGGGGTTTT" + TRIM_ADAPTER_R2),                                 # 2 full adapters
        ("AAAACCCC" + TRIM_ADAPTER_R1,
         "TTTTGGGG" + TRIM_ADAPTER_R2),                                     # 3 too short after trim
        ("AAAACCCCGGGGTT" + TRIM_ADAPTER_R1 + "GGGG",
         "TTTTGGGGCCCCAA" + TRIM_ADAPTER_R2 + "CCCC"),                      # 4 internal adapter
        ("AAAACCCCGGGG" + TRIM_ADAPTER_R1,
         "TTTTGGGGCCCCAAAATT"),                                             # 5 R1 only
        ("AAAACCCCAAAAAA" + TRIM_ADAPTER_R1[:6],
         "TTTTGGGGTTTTTT" + TRIM_ADAPTER_R2[:6]),                           # 6 partial ligation
        ("NNNNAAAACCCC" + TRIM_ADAPTER_R1,
         "NNNNTTTTGGGG" + TRIM_ADAPTER_R2),                                 # 7 undetermined
    ]
    for idx, (name, (i7, i5)) in enumerate(zip(names, index_pairs), start=1):
        base = name.split(" ")[0], name.split(" ", 1)[1]
        r1seq, r2seq = reads[idx - 1]
        r1.append(rec(f"{base[0]}/1 {base[1]}", r1seq))
        r2.append(rec(f"{base[0]}/2 {base[1]}", r2seq))
        i1.append(rec(name, i7))
        i2.append(rec(name, i5))

    _write_lane(os.path.join(out_dir, "R1.fastq.gz"), r1)
    _write_lane(os.path.join(out_dir, "R2.fastq.gz"), r2)
    _write_lane(os.path.join(out_dir, "I1.fastq.gz"), i1)
    _write_lane(os.path.join(out_dir, "I2.fastq.gz"), i2)
    _write_sheet(out_dir)

    with open(os.path.join(out_dir, "TRIMMING_README.txt"), "w", encoding="utf-8") as fh:
        fh.write(
            "Adapter trimming demo parameters (cutadapt -j 2, default matching):\n"
            f"  adapter_r1 = {TRIM_ADAPTER_R1}\n"
            f"  adapter_r2 = {TRIM_ADAPTER_R2}\n"
            "  min_length = 10\n"
            "Expected: SEQ-0003 pair is filtered (8 bp insert on both ends);\n"
            "all other pairs survive with their 3' adapters removed.\n"
        )


def main(root: str) -> None:
    if os.path.isdir(root):
        shutil.rmtree(root)
    normal_dir = os.path.join(root, "normal")
    ambiguous_dir = os.path.join(root, "ambiguous")
    corrupt_dir = os.path.join(root, "corrupt")
    trimming_dir = os.path.join(root, "trimming")
    for d in (normal_dir, ambiguous_dir, corrupt_dir, trimming_dir):
        os.makedirs(d, exist_ok=True)
    normal_dataset(normal_dir)
    ambiguous_dataset(ambiguous_dir)
    corrupt_dataset(corrupt_dir, normal_dir)
    trimming_dataset(trimming_dir)
    print(f"sample datasets written under {root}")


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ROOT
    main(os.path.abspath(root))
