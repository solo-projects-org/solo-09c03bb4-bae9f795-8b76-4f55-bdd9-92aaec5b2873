"""Tests for the cutadapt adapter-trimming feature.

Covers validation, the real cutadapt end-to-end path, resubmit/supersede
semantics, cancellation, package replacement, download immutability and the
"original demux outputs stay untouched" guarantee.
"""

from __future__ import annotations

import gzip
import io
import os
import subprocess
import threading
import time
import zipfile
from pathlib import Path

import pytest

from app.errors import TrimmingError
from app.main import app
from app.trim_jobs import (
    CANCELLED,
    COMPLETED,
    FAILED,
    QUEUED,
    RUNNING,
    SUPERSEDED,
    TrimManager,
)
from app.trimming import TrimParams, parse_report, validate_params

SAMPLES_DIR = Path(__file__).resolve().parent.parent / "samples"
ADAPTER_R1 = "AGATCGGAAGAG"
ADAPTER_R2 = "CTGTCTCTTATA"


# --------------------------------------------------------------- helpers


def _post_job(client, dataset: str):
    base = SAMPLES_DIR / dataset
    payload = []
    for key in ("R1", "R2", "I1", "I2"):
        payload.append(
            (key, (f"{key}.fastq.gz", (base / f"{key}.fastq.gz").read_bytes(),
                   "application/gzip"))
        )
    payload.append(
        ("sample_sheet", ("samples.csv", (base / "samples.csv").read_bytes(),
                          "text/csv"))
    )
    return client.post("/api/jobs", files=payload)


def _wait_attempt_terminal(client, attempt_id: str, timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        snap = client.get(f"/api/trimming/{attempt_id}").json()
        if snap["status"] in (COMPLETED, FAILED, CANCELLED, SUPERSEDED):
            return snap
        time.sleep(0.05)
    raise AssertionError("trimming attempt did not finish in time")


def _completed_trimming_batch(client) -> str:
    resp = _post_job(client, "trimming")
    assert resp.status_code == 201, resp.text
    job_id = resp.json()["id"]
    deadline = time.time() + 20
    while time.time() < deadline:
        snap = client.get(f"/api/jobs/{job_id}").json()
        if snap["status"] in ("completed", "failed", "cancelled"):
            assert snap["status"] == "completed", snap
            return job_id
        time.sleep(0.05)
    raise AssertionError("demux batch did not complete")


def _submit_trim(client, job_id: str, sample_id: str,
                 adapter_r1=ADAPTER_R1, adapter_r2=ADAPTER_R2, min_length=10):
    return client.post(
        "/api/trimming",
        data={
            "job_id": job_id,
            "sample_id": sample_id,
            "adapter_r1": adapter_r1,
            "adapter_r2": adapter_r2,
            "min_length": min_length,
        },
    )


# ------------------------------------------------------------ validation


def test_validate_params_accepts_and_normalises():
    p = validate_params(" agatcgg ", "CTGTCTCTTATA", 15)
    assert p.adapter_r1 == "AGATCGG"
    assert p.adapter_r2 == "CTGTCTCTTATA"
    assert p.min_length == 15


@pytest.mark.parametrize(
    "a1,a2,m",
    [
        ("", ADAPTER_R2, 10),
        (ADAPTER_R1, "BAD SEQ!", 10),
        (ADAPTER_R1, ADAPTER_R2, 0),
        (ADAPTER_R1, ADAPTER_R2, -3),
        (ADAPTER_R1, "AAA", 10 ** 7),
    ],
)
def test_validate_params_rejects_bad_input(a1, a2, m):
    with pytest.raises(TrimmingError):
        validate_params(a1, a2, m)


def test_parse_report_headline_counters():
    report = (
        "=== Summary ===\n"
        "Total read pairs processed:                  1,200\n"
        "Read 1 with adapter:                         900 (75.0%)\n"
        "Read 2 with adapter:                         850 (70.8%)\n"
        "== Read fate breakdown ==\n"
        "Pairs that were too short:                   100 (8.3%)\n"
        "Pairs written (passing filters):             1,100 (91.7%)\n"
    )
    stats = parse_report(report)
    assert stats == {
        "pairs_processed": 1200,
        "pairs_written": 1100,
        "pairs_filtered": 100,
        "r1_with_adapter": 900,
        "r2_with_adapter": 850,
    }


# ----------------------------------------------------------- API: targets


def test_targets_lists_only_non_empty_samples_of_completed_batches(client):
    job_id = _completed_trimming_batch(client)
    resp = client.get("/api/trimming/targets")
    assert resp.status_code == 200
    targets = resp.json()["targets"]
    keys = {(t["job_id"], t["sample_id"]) for t in targets}
    # trimming dataset: A=2, B=2, C=2 pairs; undetermined=1 -> all non-empty
    assert (job_id, "sampleA") in keys
    assert (job_id, "sampleB") in keys
    assert (job_id, "sampleC") in keys
    assert (job_id, "undetermined") in keys
    for target in targets:
        assert target["read_count"] > 0
        assert target["state"] is None  # nothing submitted yet


def test_submit_unknown_job_or_empty_sample_is_404(client):
    job_id = _completed_trimming_batch(client)
    resp = _submit_trim(client, "deadbeef", "sampleA")
    assert resp.status_code == 404
    resp = _submit_trim(client, job_id, "ghost")
    assert resp.status_code == 404


def test_submit_bad_params_is_400(client):
    job_id = _completed_trimming_batch(client)
    resp = _submit_trim(client, job_id, "sampleA", adapter_r1="nope!!")
    assert resp.status_code == 400
    assert "adapter" in resp.json()["detail"].lower()


# ----------------------------------------------------- real cutadapt path


def test_cutadapt_trims_adapters_and_filters_short_pairs(client):
    job_id = _completed_trimming_batch(client)
    resp = _submit_trim(client, job_id, "sampleB")
    assert resp.status_code == 201, resp.text
    snap = _wait_attempt_terminal(client, resp.json()["id"])
    assert snap["status"] == COMPLETED, snap
    assert snap["published"] is True
    assert snap["download_available"] is True
    # sampleB had 2 input pairs; SEQ-0003 is below min_length after trimming.
    assert snap["pairs_in"] == 2
    assert snap["pairs_out"] == 1

    dl = client.get(f"/api/trimming/{snap['id']}/download")
    assert dl.status_code == 200
    zf = zipfile.ZipFile(io.BytesIO(dl.content))
    names = set(zf.namelist())
    assert names == {
        "trimmed_R1.fastq.gz",
        "trimmed_R2.fastq.gz",
        "trim_params.txt",
        "cutadapt_report.txt",
    }
    assert zf.testzip() is None

    with gzip.open(io.BytesIO(zf.read("trimmed_R1.fastq.gz")), "rb") as fh:
        r1 = fh.read().decode()
    with gzip.open(io.BytesIO(zf.read("trimmed_R2.fastq.gz")), "rb") as fh:
        r2 = fh.read().decode()
    # The surviving pair is SEQ-0004: internal adapter means everything from
    # the adapter onward (adapter + trailing bases) is removed.
    assert "@SEQ-0004" in r1
    assert "AGATCGGAAGAG" not in r1
    assert "CTGTCTCTTATA" not in r2
    # The short pair SEQ-0003 was dropped from both ends.
    assert "SEQ-0003" not in r1 and "SEQ-0003" not in r2
    params_txt = zf.read("trim_params.txt").decode()
    assert f"adapter_r1={ADAPTER_R1}" in params_txt
    assert "pair_filter=any" in params_txt
    assert "cutadapt_threads=2" in params_txt


def test_undetermined_sample_is_trimmable_too(client):
    job_id = _completed_trimming_batch(client)
    resp = _submit_trim(client, job_id, "undetermined")
    assert resp.status_code == 201
    snap = _wait_attempt_terminal(client, resp.json()["id"])
    assert snap["status"] == COMPLETED, snap
    assert snap["pairs_out"] == 1


def test_demux_originals_are_never_modified(client):
    job_id = _completed_trimming_batch(client)
    from app.main import manager

    job = manager.get(job_id)
    r1_path = job.work_dir / "output" / "sampleB_R1.fastq.gz"
    before = r1_path.read_bytes()
    resp = _submit_trim(client, job_id, "sampleB")
    _wait_attempt_terminal(client, resp.json()["id"])
    assert r1_path.read_bytes() == before
    # Original demux download still works and still contains the adapter.
    dl = client.get(f"/api/jobs/{job_id}/download")
    assert dl.status_code == 200
    zf = zipfile.ZipFile(io.BytesIO(dl.content))
    with gzip.open(io.BytesIO(zf.read("sampleB_R1.fastq.gz")), "rb") as fh:
        assert ADAPTER_R1.encode() in fh.read()


# ------------------------------------------------------- supersede/cancel


def test_resubmit_supersedes_running_attempt_and_old_never_publishes(
    client, monkeypatch
):
    import app.trim_jobs as trim_jobs_module

    job_id = _completed_trimming_batch(client)
    entered = threading.Event()
    release = threading.Event()
    real_run = trim_jobs_module.run_trimming

    def gated_run(source_r1=None, source_r2=None, work_dir=None, params=None,
                  is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        entered.set()
        # Stay here until either cancelled (superseded) or released.
        while not release.wait(timeout=0.02):
            if is_cancelled():
                from app.errors import DemuxError

                raise DemuxError("trimming cancelled by user")
        return real_run(source_r1, source_r2, work_dir, params, is_cancelled,
                        launch_subprocess)

    monkeypatch.setattr(trim_jobs_module, "run_trimming", gated_run)
    try:
        first = _submit_trim(client, job_id, "sampleA", min_length=10)
        assert entered.wait(timeout=10)
        # Resubmit with different parameters: must be accepted immediately.
        second = _submit_trim(client, job_id, "sampleA", min_length=12)
        assert second.status_code == 201
        first_snap = _wait_attempt_terminal(client, first.json()["id"])
        assert first_snap["status"] == SUPERSEDED
        assert first_snap["published"] is False
        assert first_snap["download_available"] is False

        # Old attempt's temporary directory must be gone before the new one
        # starts publishing; release the gate and let the last attempt run.
        release.set()
        second_snap = _wait_attempt_terminal(client, second.json()["id"])
        assert second_snap["status"] == COMPLETED
        assert second_snap["params"]["min_length"] == 12
        assert second_snap["published"] is True
    finally:
        release.set()
        monkeypatch.undo()


def test_rapid_resubmits_retire_displaced_queued_attempt(client, monkeypatch):
    """A queued attempt replaced while the worker is busy must not be orphaned."""
    import app.trim_jobs as trim_jobs_module
    from app.main import TRIM_ROOT

    job_id = _completed_trimming_batch(client)
    entered = threading.Event()
    release = threading.Event()
    real_run = trim_jobs_module.run_trimming

    def gated_run(source_r1=None, source_r2=None, work_dir=None, params=None,
                  is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        entered.set()
        while not release.wait(timeout=0.02):
            if is_cancelled():
                from app.errors import DemuxError

                raise DemuxError("trimming cancelled by user")
        return real_run(source_r1, source_r2, work_dir, params, is_cancelled,
                        launch_subprocess)

    monkeypatch.setattr(trim_jobs_module, "run_trimming", gated_run)
    try:
        first = _submit_trim(client, job_id, "sampleA", min_length=10)
        assert entered.wait(timeout=10)
        # The worker is busy: the second submission only queues.
        second = _submit_trim(client, job_id, "sampleA", min_length=11)
        assert second.json()["status"] in (QUEUED, RUNNING)
        second_id = second.json()["id"]
        second_dir = TRIM_ROOT / trim_jobs_module.TrimManager.target_key(
            job_id, "sampleA"
        ).replace(":", "_") / second_id
        # A third submission immediately displaces the queued second one.
        third = _submit_trim(client, job_id, "sampleA", min_length=12)
        third_id = third.json()["id"]

        # The displaced queued attempt must already be terminal (not stuck in
        # queued/running), never published, and its temp directory removed.
        second_snap = client.get(f"/api/trimming/{second_id}").json()
        assert second_snap["status"] == SUPERSEDED
        assert second_snap["published"] is False
        deadline = time.time() + 10
        while time.time() < deadline and second_dir.exists():
            time.sleep(0.05)
        assert not second_dir.exists()
        # Cancelling an already-terminal attempt returns immediately.
        cancel = client.post(f"/api/trimming/{second_id}/cancel")
        assert cancel.status_code == 200
        assert cancel.json()["status"] == SUPERSEDED

        release.set()
        first_snap = _wait_attempt_terminal(client, first.json()["id"])
        third_snap = _wait_attempt_terminal(client, third_id)
        assert first_snap["status"] == SUPERSEDED
        assert third_snap["status"] == COMPLETED
        assert third_snap["published"] is True
        assert third_snap["params"]["min_length"] == 12
    finally:
        release.set()
        monkeypatch.undo()


def test_cancel_waits_for_process_exit_and_keeps_old_package(
    client, monkeypatch
):
    import app.trim_jobs as trim_jobs_module

    job_id = _completed_trimming_batch(client)
    # Establish a first, genuinely published package.
    first = _submit_trim(client, job_id, "sampleC", min_length=10)
    first_snap = _wait_attempt_terminal(client, first.json()["id"])
    assert first_snap["status"] == COMPLETED

    entered = threading.Event()
    release = threading.Event()
    real_run = trim_jobs_module.run_trimming

    def gated_run(source_r1=None, source_r2=None, work_dir=None, params=None,
                  is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        entered.set()
        while not release.wait(timeout=0.02):
            if is_cancelled():
                from app.errors import DemuxError

                raise DemuxError("trimming cancelled by user")
        return real_run(source_r1, source_r2, work_dir, params, is_cancelled,
                        launch_subprocess)

    monkeypatch.setattr(trim_jobs_module, "run_trimming", gated_run)
    try:
        second = _submit_trim(client, job_id, "sampleC", min_length=8)
        assert entered.wait(timeout=10)
        cancel_resp = client.post(
            f"/api/trimming/{second.json()['id']}/cancel"
        )
        assert cancel_resp.status_code == 200
        snap = cancel_resp.json()
        assert snap["status"] == CANCELLED
        # The cancel endpoint returned only after the worker stopped.
        assert snap["finished_at"] is not None
        # Old package remains the published one and stays downloadable.
        assert client.get(
            f"/api/trimming/{first_snap['id']}/download"
        ).status_code == 200
        # Cancelled attempt cannot be downloaded.
        assert client.get(
            f"/api/trimming/{second.json()['id']}/download"
        ).status_code == 404
    finally:
        release.set()
        monkeypatch.undo()


def test_failed_attempt_keeps_old_package_and_is_marked(client, monkeypatch):
    import app.trim_jobs as trim_jobs_module

    job_id = _completed_trimming_batch(client)
    first = _submit_trim(client, job_id, "sampleA", min_length=10)
    first_snap = _wait_attempt_terminal(client, first.json()["id"])
    assert first_snap["status"] == COMPLETED

    def boom(source_r1=None, source_r2=None, work_dir=None, params=None,
             is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        raise TrimmingError("simulated cutadapt failure")

    monkeypatch.setattr(trim_jobs_module, "run_trimming", boom)
    try:
        second = _submit_trim(client, job_id, "sampleA", min_length=11)
        snap = _wait_attempt_terminal(client, second.json()["id"])
        assert snap["status"] == FAILED
        assert "simulated cutadapt failure" in snap["error"]
        assert snap["download_available"] is False
        # Old package still published.
        assert client.get(
            f"/api/trimming/{first_snap['id']}/download"
        ).status_code == 200
        listed = client.get("/api/trimming").json()["attempts"]
        failed = next(a for a in listed if a["id"] == second.json()["id"])
        assert failed["params"]["min_length"] == 11
        # The failed attempt must point back at the still-active old package.
        assert failed["active_published"] is not None
        assert failed["active_published"]["attempt_id"] == first_snap["id"]
        assert failed["active_published"]["params"]["min_length"] == 10
        first_listed = next(a for a in listed if a["id"] == first_snap["id"])
        assert first_listed["active_published"] is None
    finally:
        monkeypatch.undo()


def test_cancel_does_not_block_other_samples(client, monkeypatch):
    import app.trim_jobs as trim_jobs_module

    job_id = _completed_trimming_batch(client)
    gate_a = threading.Event()
    entered_a = threading.Event()
    real_run = trim_jobs_module.run_trimming

    def gate_sample_a(source_r1=None, source_r2=None, work_dir=None, params=None,
                      is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        # Only block sampleA; sampleC runs on its own worker thread.
        if params.adapter_r1 == ADAPTER_R1 and "sampleA" in work_dir:
            entered_a.set()
            gate_a.wait(timeout=15)
        return real_run(source_r1, source_r2, work_dir, params, is_cancelled,
                        launch_subprocess)

    monkeypatch.setattr(trim_jobs_module, "run_trimming", gate_sample_a)
    try:
        a = _submit_trim(client, job_id, "sampleA")
        assert entered_a.wait(timeout=10)
        # sampleC is a different target: its worker must run straight through.
        c = _submit_trim(client, job_id, "sampleC")
        c_snap = _wait_attempt_terminal(client, c.json()["id"])
        assert c_snap["status"] == COMPLETED

        cancel = client.post(f"/api/trimming/{a.json()['id']}/cancel")
        assert cancel.status_code == 200
        assert cancel.json()["status"] == CANCELLED
    finally:
        gate_a.set()
        monkeypatch.undo()


# ------------------------------------------------------------- downloads


def test_superseded_and_unknown_attempts_block_download(client):
    job_id = _completed_trimming_batch(client)
    assert client.get("/api/trimming/deadbeef/download").status_code == 404
    # A queued/running attempt has no ZIP yet.
    resp = _submit_trim(client, job_id, "sampleB")
    attempt_id = resp.json()["id"]
    # Immediately after creation it may already be running/queued; either way
    # no download until completion.
    snap = client.get(f"/api/trimming/{attempt_id}").json()
    if snap["status"] in (QUEUED, RUNNING):
        assert client.get(
            f"/api/trimming/{attempt_id}/download"
        ).status_code == 404
    _wait_attempt_terminal(client, attempt_id)


def test_inflight_download_receives_one_complete_version(client, monkeypatch):
    """A package replacement must not truncate a download that already started."""
    job_id = _completed_trimming_batch(client)
    first = _submit_trim(client, job_id, "sampleA", min_length=10)
    first_snap = _wait_attempt_terminal(client, first.json()["id"])
    assert first_snap["status"] == COMPLETED

    # Open the first package the same way the download endpoint does.
    from app.main import trim_manager

    fh, _name = trim_manager.open_published_zip(first_snap["id"])
    attempt_id = first_snap["id"]
    assert fh.read(4) == b"PK\x03\x04"
    fh.seek(0)

    # Replace the published package with new parameters.
    second = _submit_trim(client, job_id, "sampleA", min_length=13)
    second_snap = _wait_attempt_terminal(client, second.json()["id"])
    assert second_snap["status"] == COMPLETED

    # The old fd keeps serving the old, now-unlinked file to completion.
    data = fh.read()
    fh.close()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        assert zf.testzip() is None
        assert "trimmed_R1.fastq.gz" in zf.namelist()
    # The old attempt is no longer downloadable through the API.
    assert client.get(f"/api/trimming/{attempt_id}/download").status_code == 404
    assert client.get(
        f"/api/trimming/{second_snap['id']}/download"
    ).status_code == 200


def test_temp_directories_are_cleaned_after_supersede(client, monkeypatch):
    import app.trim_jobs as trim_jobs_module
    from app.main import TRIM_ROOT

    job_id = _completed_trimming_batch(client)
    entered = threading.Event()
    release = threading.Event()
    real_run = trim_jobs_module.run_trimming

    def gated_run(source_r1=None, source_r2=None, work_dir=None, params=None,
                  is_cancelled=None, launch_subprocess=subprocess.Popen, **_kw):
        entered.set()
        while not release.wait(timeout=0.02):
            if is_cancelled():
                from app.errors import DemuxError

                raise DemuxError("trimming cancelled by user")
        return real_run(source_r1, source_r2, work_dir, params, is_cancelled,
                        launch_subprocess)

    monkeypatch.setattr(trim_jobs_module, "run_trimming", gated_run)
    try:
        first = _submit_trim(client, job_id, "sampleB")
        assert entered.wait(timeout=10)
        first_id = first.json()["id"]
        first_dir = TRIM_ROOT / trim_jobs_module.TrimManager.target_key(
            job_id, "sampleB"
        ).replace(":", "_") / first_id
        assert first_dir.exists()
        second = _submit_trim(client, job_id, "sampleB", min_length=14)
        _wait_attempt_terminal(client, first_id)
        # Old attempt fully cleaned before/at supersede completion.
        deadline = time.time() + 10
        while time.time() < deadline and first_dir.exists():
            time.sleep(0.05)
        assert not first_dir.exists()
        release.set()
        _wait_attempt_terminal(client, second.json()["id"])
    finally:
        release.set()
        monkeypatch.undo()


# ----------------------------------------------------- unit: process group


def test_run_trimming_uses_cutadapt_j2_pair_filter_and_3p_adapters(tmp_path):
    import subprocess

    from app import trimming as trimming_module

    captured: dict = {}

    class FakePopen:
        def __init__(self, argv, **kwargs):
            self.argv = argv
            captured["argv"] = argv
            captured["kwargs"] = kwargs
            # cutadapt's report is written through the stdout file handle.
            report_fh = kwargs["stdout"]
            report_fh.write(
                b"Total read pairs processed: 1\n"
                b"Pairs written (passing filters): 1\n"
            )
            self.returncode = 0
            self.pid = 4242

        def communicate(self, timeout=None):
            # Simulate cutadapt writing the outputs named in the argv.
            out_r1 = self.argv[self.argv.index("-o") + 1]
            out_r2 = self.argv[self.argv.index("-p") + 1]
            for name in (out_r1, out_r2):
                Path(name).parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(name, "wb") as fh:
                    fh.write(b"@r1\nACGTACGT\n+\nIIIIIIII\n")
            return None, b""

        def poll(self):
            return 0

    params = validate_params(ADAPTER_R1, ADAPTER_R2, 10)
    source = tmp_path / "src"
    source.mkdir()
    r1 = source / "s_R1.fastq.gz"
    r2 = source / "s_R2.fastq.gz"
    for p in (r1, r2):
        with gzip.open(p, "wb") as fh:
            fh.write(b"@r1\nACGTACGT\n+\nIIIIIIII\n")
    result = trimming_module.run_trimming(
        str(r1), str(r2), str(tmp_path / "work"), params,
        lambda: False, launch_subprocess=FakePopen,
    )
    argv = captured["argv"]
    assert "-j" in argv and argv[argv.index("-j") + 1] == "2"
    assert argv[argv.index("-a") + 1] == ADAPTER_R1
    assert argv[argv.index("-A") + 1] == ADAPTER_R2
    assert "--pair-filter=any" in argv
    assert "-m" in argv and argv[argv.index("-m") + 1] == "10"
    assert "--report=full" in argv
    assert captured["kwargs"]["start_new_session"] is True
    assert hasattr(captured["kwargs"]["stdout"], "write")
    assert Path(result.zip_path).is_file()


def test_run_trimming_nonzero_exit_becomes_trimming_error(tmp_path):
    from app import trimming as trimming_module

    class FailingPopen:
        def __init__(self, argv, **kwargs):
            self.returncode = 1
            self.pid = 1

        def communicate(self, timeout=None):
            return b"", b"cutadapt: bad adapter"

        def poll(self):
            return 1

    params = validate_params(ADAPTER_R1, ADAPTER_R2, 10)
    source = tmp_path / "src"
    source.mkdir()
    r1 = source / "s_R1.fastq.gz"
    r2 = source / "s_R2.fastq.gz"
    for p in (r1, r2):
        with gzip.open(p, "wb") as fh:
            fh.write(b"@r1\nACGT\n+\nIIII\n")
    with pytest.raises(TrimmingError, match="cutadapt exited"):
        trimming_module.run_trimming(
            str(r1), str(r2), str(tmp_path / "work"), params,
            lambda: False, launch_subprocess=FailingPopen,
        )


def test_run_trimming_cancellation_kills_process_group(tmp_path):
    import signal

    from app import trimming as trimming_module

    class SlowPopen:
        def __init__(self, argv, **kwargs):
            self.argv = argv
            self.pid = os.getpid()
            self.returncode = None
            self.signals: list[int] = []

        def communicate(self, timeout=None):
            # Block until the termination helper marks the process dead.
            if self.returncode is None:
                raise subprocess.TimeoutExpired(self.argv, timeout)
            return None, b""

        def poll(self):
            return self.returncode

        def terminate(self):
            self.signals.append(signal.SIGTERM)
            self.returncode = -15

        def kill(self):
            self.signals.append(signal.SIGKILL)
            self.returncode = -9

    proc_holder: dict[str, SlowPopen] = {}
    seen_signals: list[int] = []
    real_killpg = trimming_module.os.killpg

    def fake_killpg(pgid, sig):
        seen_signals.append(sig)
        proc = proc_holder.get("proc")
        if proc is not None and sig == signal.SIGTERM:
            proc.terminate()

    def launch(argv, **kwargs):
        proc = SlowPopen(argv, **kwargs)
        proc_holder["proc"] = proc
        return proc

    trimming_module.os.killpg = fake_killpg
    # Flip to cancelled only after the (fake) process has been launched.
    cancel_event = threading.Event()

    def is_cancelled() -> bool:
        return "proc" in proc_holder and cancel_event.is_set()

    def flip_after_launch():
        deadline = time.time() + 5
        while "proc" not in proc_holder and time.time() < deadline:
            time.sleep(0.01)
        cancel_event.set()

    threading.Thread(target=flip_after_launch, daemon=True).start()
    try:
        params = validate_params(ADAPTER_R1, ADAPTER_R2, 10)
        source = tmp_path / "src"
        source.mkdir()
        r1 = source / "s_R1.fastq.gz"
        r2 = source / "s_R2.fastq.gz"
        for p in (r1, r2):
            with gzip.open(p, "wb") as fh:
                fh.write(b"@r1\nACGT\n+\nIIII\n")
        from app.errors import DemuxError

        with pytest.raises(DemuxError, match="cancelled"):
            trimming_module.run_trimming(
                str(r1), str(r2), str(tmp_path / "work"), params,
                is_cancelled,
                launch_subprocess=launch,
            )
        # The whole -j 2 process group must have been asked to stop.
        assert signal.SIGTERM in seen_signals
        assert proc_holder["proc"].returncode == -15
    finally:
        trimming_module.os.killpg = real_killpg
