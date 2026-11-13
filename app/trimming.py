"""Adapter-trimming pipeline built on top of a published demux batch.

The demultiplex output (``<sample>_R1.fastq.gz`` / ``<sample>_R2.fastq.gz``)
is *never* modified: cutadapt reads those files and writes a brand-new pair
into the attempt's private work directory.  cutadapt is invoked with its
default 3' adapter matching rules, two worker threads (``-j 2``) and the
"any" pair filter so that a pair is discarded whenever **either** end ends
up shorter than the minimum length.

Only a fully finished, verified pair plus a complete ZIP is handed back to
the caller; every failure or cancellation raises and leaves publishing to
the manager layer.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import zipfile
from dataclasses import dataclass
from typing import Callable

from .errors import DemuxError, TrimmingError

# cutadapt's own default error tolerance / overlap; we pass nothing that
# changes matching behaviour, these constants only document the defaults
# the task relies on.
CUTADAPT_THREADS = 2
MIN_ADAPTER_LEN = 1
MAX_ADAPTER_LEN = 100
MIN_READ_LENGTH = 1
MAX_READ_LENGTH = 1_000_000
_ALLOWED_ADAPTER_CHARS = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_OUTPUT_SUFFIXES = ("R1", "R2")


@dataclass(frozen=True, slots=True)
class TrimParams:
    """User-submitted trimming parameters for one sample."""

    adapter_r1: str
    adapter_r2: str
    min_length: int

    def as_dict(self) -> dict[str, object]:
        return {
            "adapter_r1": self.adapter_r1,
            "adapter_r2": self.adapter_r2,
            "min_length": self.min_length,
        }

    def command_line(self, r1_in: str, r2_in: str, r1_out: str, r2_out: str,
                     report_path: str) -> list[str]:
        """Return the exact cutadapt argv used for this attempt.

        cutadapt 5.x writes its report to stdout (``--report`` only selects
        full/minimal), so the caller redirects stdout into ``report_path``.
        """
        return [
            "cutadapt",
            "-j",
            str(CUTADAPT_THREADS),
            # 3' adapters; default anchored-overlap / error-rate rules apply.
            "-a",
            self.adapter_r1,
            "-A",
            self.adapter_r2,
            # Drop the whole pair when either end is shorter than the limit.
            # "any" is cutadapt's paired-end default; state it explicitly so
            # the requirement is self-documenting.
            "--pair-filter=any",
            "-m",
            str(self.min_length),
            "-o",
            r1_out,
            "-p",
            r2_out,
            "--report=full",
            r1_in,
            r2_in,
        ]


def validate_params(adapter_r1: str, adapter_r2: str, min_length: int) -> TrimParams:
    """Validate the web form values and return a frozen parameter object."""

    adapters: dict[str, str] = {}
    for label, raw in (("adapter_r1", adapter_r1), ("adapter_r2", adapter_r2)):
        if not isinstance(raw, str):
            raise TrimmingError(f"{label} must be a string")
        seq = raw.strip().upper()
        if not (MIN_ADAPTER_LEN <= len(seq) <= MAX_ADAPTER_LEN):
            raise TrimmingError(
                f"{label} length must be between {MIN_ADAPTER_LEN} and "
                f"{MAX_ADAPTER_LEN} bases"
            )
        bad = sorted(set(seq) - _ALLOWED_ADAPTER_CHARS)
        if bad:
            raise TrimmingError(
                f"{label} may only contain IUPAC nucleotide codes "
                f"(A/C/G/T/N/...), got unexpected characters: {''.join(bad)!r}"
            )
        adapters[label] = seq

    if isinstance(min_length, bool) or not isinstance(min_length, int):
        raise TrimmingError("min_length must be an integer")
    if not (MIN_READ_LENGTH <= min_length <= MAX_READ_LENGTH):
        raise TrimmingError(
            f"min_length must be between {MIN_READ_LENGTH} and {MAX_READ_LENGTH}"
        )

    return TrimParams(
        adapter_r1=adapters["adapter_r1"],
        adapter_r2=adapters["adapter_r2"],
        min_length=min_length,
    )


def _resolve_cutadapt() -> list[str]:
    """Locate the cutadapt executable (PATH, then ``python -m cutadapt``)."""

    exe = shutil.which("cutadapt")
    if exe:
        return [exe]
    import sys

    return [sys.executable, "-m", "cutadapt"]


def _count_fastq_records(path: str) -> int:
    """Count FASTQ records in a gzip file, validating the 4-line structure."""

    count = 0
    with gzip.open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not header.startswith(b"@") or not plus.startswith(b"+") or not seq or not qual:
                raise TrimmingError(f"cutadapt output {os.path.basename(path)} is not valid FASTQ")
            if len(seq.rstrip(b"\n")) != len(qual.rstrip(b"\n")):
                raise TrimmingError(
                    f"cutadapt output {os.path.basename(path)} has mismatched "
                    "sequence/quality lengths"
                )
            count += 1
    return count


# cutadapt's summary lines look like:
#   Total read pairs processed:                  7
#   Pairs written (passing filters):             5 (71.4%)
#   Pairs that were too short:                   2 (28.6%)
#   Read 1 with adapter:                         6 (85.7%)
_NUM_RE = re.compile(r"([\d,]+)")


def _first_int(pattern: str, text: str) -> int | None:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return None
    # Prefer an explicit capture group (so "Read 1 with adapter: 900"
    # returns 900 rather than the digit inside "Read 1").
    token = match.group(1) if match.groups() else _NUM_RE.search(match.group(0)).group(1)
    return int(token.replace(",", ""))


def parse_report(report: str) -> dict[str, int | None]:
    """Pull the headline counters out of a cutadapt text report."""

    return {
        "pairs_processed": _first_int(
            r"(?:Total )?read pairs processed:\s*([\d,]+)", report
        ),
        "pairs_written": _first_int(
            r"Pairs written[^:\n]*:\s*([\d,]+)", report
        ),
        "pairs_filtered": _first_int(
            r"Pairs (?:that were too short|filtered)[^:\n]*:\s*([\d,]+)", report
        ),
        "r1_with_adapter": _first_int(
            r"Read 1 with adapter:\s*([\d,]+)", report
        ),
        "r2_with_adapter": _first_int(
            r"Read 2 with adapter:\s*([\d,]+)", report
        ),
    }


@dataclass(slots=True)
class TrimResult:
    zip_path: str
    output_r1: str
    output_r2: str
    report_path: str
    stats: dict[str, int | None]
    pairs_in: int
    pairs_out: int


def run_trimming(
    source_r1: str,
    source_r2: str,
    work_dir: str,
    params: TrimParams,
    is_cancelled: Callable[[], bool],
    launch_subprocess: Callable[..., subprocess.Popen] = subprocess.Popen,
) -> TrimResult:
    """Run one cutadapt attempt and verify + package its outputs.

    ``launch_subprocess`` is injectable so tests can run without a real
    cutadapt binary.  Raises :class:`TrimmingError` on any failure and
    :class:`DemuxError` on cancellation.
    """

    for label, path in (("source R1", source_r1), ("source R2", source_r2)):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise TrimmingError(f"{label} file is missing or empty: {path}")

    os.makedirs(work_dir, exist_ok=True)
    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    r1_out = os.path.join(output_dir, "trimmed_R1.fastq.gz")
    r2_out = os.path.join(output_dir, "trimmed_R2.fastq.gz")
    report_path = os.path.join(work_dir, "cutadapt_report.txt")
    argv = params.command_line(source_r1, source_r2, r1_out, r2_out, report_path)
    if argv[0] == "cutadapt":
        argv[0:1] = _resolve_cutadapt()

    if is_cancelled():
        raise DemuxError("trimming cancelled by user")

    try:
        # Own process group so -j 2 child workers can all be signalled at once.
        # cutadapt writes its human-readable report to stdout; capture it in
        # the report file while stderr stays in memory for error messages.
        report_fh = open(report_path, "wb")  # noqa: SIM115 - closed below
        try:
            process = launch_subprocess(
                argv,
                stdout=report_fh,
                stderr=subprocess.PIPE,
                start_new_session=True,
            )
        except OSError as exc:
            report_fh.close()
            raise TrimmingError(f"failed to start cutadapt: {exc}") from exc
    except OSError as exc:
        raise TrimmingError(f"failed to start cutadapt: {exc}") from exc

    try:
        while True:
            try:
                _stdout, stderr = process.communicate(timeout=0.2)
                break
            except subprocess.TimeoutExpired:
                if is_cancelled():
                    _terminate_process_group(process)
                    try:
                        process.communicate(timeout=30)
                    except subprocess.TimeoutExpired:
                        _kill_process_group(process)
                        process.communicate()
                    raise DemuxError("trimming cancelled by user")
    except DemuxError:
        raise
    finally:
        # Make sure no cutadapt process ever outlives an attempt.
        if process.poll() is None:
            _terminate_process_group(process)
        report_fh.close()

    if process.returncode != 0:
        detail = (stderr or b"").decode("utf-8", errors="replace").strip()
        if not detail and os.path.isfile(report_path):
            with open(report_path, "r", encoding="utf-8", errors="replace") as fh:
                detail = fh.read().strip()
        raise TrimmingError(
            f"cutadapt exited with code {process.returncode}: "
            f"{detail[-500:] or 'no error output'}"
        )

    if is_cancelled():
        raise DemuxError("trimming cancelled by user")

    # Verify the advertised outputs really exist, are gzip FASTQ and stay
    # paired (R1 and R2 must keep the same record count).
    for path in (r1_out, r2_out):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise TrimmingError(
                f"cutadapt reported success but produced no output: {os.path.basename(path)}"
            )
    pairs_r1 = _count_fastq_records(r1_out)
    pairs_r2 = _count_fastq_records(r2_out)
    if pairs_r1 != pairs_r2:
        raise TrimmingError(
            "trimmed R1/R2 record counts differ "
            f"({pairs_r1} vs {pairs_r2}); paired output is inconsistent"
        )

    report = ""
    if os.path.isfile(report_path):
        with open(report_path, "r", encoding="utf-8", errors="replace") as fh:
            report = fh.read()
    stats = parse_report(report)
    pairs_in = stats["pairs_processed"] if stats["pairs_processed"] is not None else -1
    pairs_out = stats["pairs_written"] if stats["pairs_written"] is not None else pairs_r1

    params_path = os.path.join(output_dir, "trim_params.txt")
    with open(params_path, "w", encoding="utf-8") as fh:
        fh.write("adapter_r1=" + params.adapter_r1 + "\n")
        fh.write("adapter_r2=" + params.adapter_r2 + "\n")
        fh.write(f"min_length={params.min_length}\n")
        fh.write("cutadapt_threads=" + str(CUTADAPT_THREADS) + "\n")
        fh.write("pair_filter=any\n")
        fh.write("command=" + " ".join(argv) + "\n")

    zip_path = os.path.join(work_dir, "trimmed.zip")
    _build_zip(zip_path, output_dir, report_path)

    # Final integrity gate: a truncated ZIP must never be published.
    with zipfile.ZipFile(zip_path) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise TrimmingError(f"generated ZIP failed integrity check: {bad}")

    return TrimResult(
        zip_path=zip_path,
        output_r1=r1_out,
        output_r2=r2_out,
        report_path=report_path,
        stats=stats,
        pairs_in=pairs_in,
        pairs_out=pairs_out,
    )


def _build_zip(zip_path: str, output_dir: str, report_path: str) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(os.path.join(output_dir, "trimmed_R1.fastq.gz"),
                 arcname="trimmed_R1.fastq.gz")
        zf.write(os.path.join(output_dir, "trimmed_R2.fastq.gz"),
                 arcname="trimmed_R2.fastq.gz")
        zf.write(os.path.join(output_dir, "trim_params.txt"),
                 arcname="trim_params.txt")
        if os.path.isfile(report_path):
            zf.write(report_path, arcname="cutadapt_report.txt")


def _terminate_process_group(process: subprocess.Popen) -> None:
    import signal

    try:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        try:
            process.terminate()
        except ProcessLookupError:
            pass


def _kill_process_group(process: subprocess.Popen) -> None:
    import signal

    try:
        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        try:
            process.kill()
        except ProcessLookupError:
            pass
