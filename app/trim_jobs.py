"""In-memory registry for adapter-trimming attempts.

Concurrency model
-----------------
Each trim *target* (one sample group inside one completed demux batch) owns a
single serial worker thread.  Submissions for the same sample are therefore
strictly ordered: a newly accepted parameter set supersedes whatever is
running or queued, the old cutadapt process group is terminated, its whole
attempt directory is removed, and only then does the worker start the last
accepted attempt.  Different samples - and the demux pipeline itself - run
on independent threads and are never blocked by a cancel on another sample.

State is memory-only; the trimming root is wiped at process startup.
"""

from __future__ import annotations

import shutil
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .errors import DemuxError, TrimmingError
from .trimming import TrimParams, run_trimming

# Attempt lifecycle states.
QUEUED = "queued"
RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"
CANCELLED = "cancelled"
SUPERSEDED = "superseded"
TERMINAL_STATES = frozenset((COMPLETED, FAILED, CANCELLED, SUPERSEDED))

# Cancel must wait for the real cutadapt process group to exit and for the
# attempt directory to be removed.  cutadapt shutdown is normally immediate;
# this is a defensive upper bound.
CANCEL_WAIT_TIMEOUT = 300.0


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


@dataclass
class TrimAttempt:
    """One accepted parameter submission for one target."""

    id: str
    job_id: str
    sample_id: str
    params: TrimParams
    work_dir: Path
    generation: int
    status: str = QUEUED
    error: str | None = None
    stats: dict[str, int | None] | None = None
    pairs_in: int | None = None
    pairs_out: int | None = None
    zip_path: Path | None = None
    created_at: datetime = field(default_factory=_utc_now)
    finished_at: datetime | None = None
    cancel_event: threading.Event = field(default_factory=threading.Event)
    user_cancelled: bool = False
    terminated: threading.Event = field(default_factory=threading.Event)

    def snapshot(self, published: bool = False) -> dict[str, Any]:
        return {
            "id": self.id,
            "job_id": self.job_id,
            "sample_id": self.sample_id,
            "status": self.status,
            "error": self.error,
            "params": self.params.as_dict(),
            "stats": dict(self.stats) if self.stats else None,
            "pairs_in": self.pairs_in,
            "pairs_out": self.pairs_out,
            "published": published,
            "download_available": published and self.zip_path is not None
            and self.zip_path.exists(),
            "created_at": self.created_at.isoformat(),
            "finished_at": self.finished_at.isoformat() if self.finished_at else None,
        }


class _TargetState:
    """Serialization point for all attempts touching one sample."""

    def __init__(self, key: str, job_id: str, sample_id: str, root: Path) -> None:
        self.key = key
        self.job_id = job_id
        self.sample_id = sample_id
        self.root = root
        self.cond = threading.Condition()
        self.pending: TrimAttempt | None = None
        self.current: TrimAttempt | None = None
        self.latest_generation = 0
        self.published: TrimAttempt | None = None
        self.worker_started = False
        self.closed = False

    def snapshot(self) -> dict[str, Any]:
        """A target-level view for the sample picker."""

        latest = self.current or self.pending
        published = self.published
        return {
            "target_id": self.key,
            "job_id": self.job_id,
            "sample_id": self.sample_id,
            "latest_attempt": latest.snapshot() if latest else None,
            "published_attempt": (
                published.snapshot(published=True) if published else None
            ),
        }


class TrimManager:
    def __init__(self, trimming_root: Path, demux_manager: Any) -> None:
        self.trimming_root = Path(trimming_root)
        self._demux = demux_manager
        self._lock = threading.Lock()
        self._targets: dict[str, _TargetState] = {}
        self._attempts: dict[str, TrimAttempt] = {}

    # ------------------------------------------------------------ startup

    def wipe(self) -> None:
        """Remove all on-disk trimming data (state is memory-only)."""

        shutil.rmtree(self.trimming_root, ignore_errors=True)
        self.trimming_root.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------- lookup

    @staticmethod
    def target_key(job_id: str, sample_id: str) -> str:
        return f"{job_id}:{sample_id}"

    def resolve_target(self, job_id: str, sample_id: str) -> tuple[Any, str, str]:
        """Validate that ``sample_id`` is a non-empty output of a completed batch.

        Returns ``(demux_job, source_r1, source_r2)`` or raises KeyError.
        """

        job = self._demux.get(job_id)
        if job is None or job.status != "completed":
            raise KeyError("completed demux batch not found")
        counts = job.final_counts or {}
        if sample_id not in counts or counts[sample_id] <= 0:
            raise KeyError("sample has no non-empty output in this batch")
        output_dir = job.work_dir / "output"
        r1 = output_dir / f"{sample_id}_R1.fastq.gz"
        r2 = output_dir / f"{sample_id}_R2.fastq.gz"
        if not (r1.is_file() and r2.is_file() and r1.stat().st_size > 0
                and r2.stat().st_size > 0):
            raise KeyError("sample output files are missing on disk")
        return job, str(r1), str(r2)

    def _get_or_create_target(self, job_id: str, sample_id: str) -> _TargetState:
        key = self.target_key(job_id, sample_id)
        with self._lock:
            target = self._targets.get(key)
            if target is None:
                target = _TargetState(
                    key, job_id, sample_id, self.trimming_root / key.replace(":", "_")
                )
                self._targets[key] = target
            return target

    # ----------------------------------------------------------- submission

    def submit(
        self, job_id: str, sample_id: str, params: TrimParams
    ) -> TrimAttempt:
        """Accept a parameter set, replacing any earlier attempt for the sample."""

        # Validate before touching any state.
        self.resolve_target(job_id, sample_id)
        target = self._get_or_create_target(job_id, sample_id)

        attempt_id = uuid.uuid4().hex
        work_dir = target.root / attempt_id
        work_dir.mkdir(parents=True, exist_ok=False)
        attempt = TrimAttempt(
            id=attempt_id,
            job_id=job_id,
            sample_id=sample_id,
            params=params,
            work_dir=work_dir,
            generation=0,  # assigned under the target condition below
        )
        with self._lock:
            self._attempts[attempt_id] = attempt

        displaced: list[TrimAttempt] = []
        with target.cond:
            target.latest_generation += 1
            attempt.generation = target.latest_generation

            # A resubmit always supersedes the in-flight / queued attempt.
            # The running attempt is signalled and retired by the worker when
            # cutadapt exits.  A *queued* attempt, though, would never be seen
            # by the worker once pending is overwritten, so retire it here and
            # now: mark it terminal, release anyone waiting in cancel(), and
            # remember its directory for cleanup outside the condition.
            old_pending = target.pending
            if old_pending is not None and old_pending is not attempt:
                old_pending.cancel_event.set()
                old_pending.status = (
                    CANCELLED if old_pending.user_cancelled else SUPERSEDED
                )
                old_pending.finished_at = _utc_now()
                old_pending.terminated.set()
                displaced.append(old_pending)
            if target.current is not None and target.current is not attempt:
                target.current.cancel_event.set()

            target.pending = attempt
            if not target.worker_started:
                target.worker_started = True
                thread = threading.Thread(
                    target=self._worker_loop,
                    args=(target,),
                    name=f"trim-{target.key}",
                    daemon=True,
                )
                thread.start()
            target.cond.notify_all()

        for old in displaced:
            # Never ran cutadapt; the directory only holds the empty attempt
            # folder created above.  Remove it without holding the condition.
            self._remove_tree(old.work_dir)
        return attempt

    # --------------------------------------------------------------- worker

    def _worker_loop(self, target: _TargetState) -> None:
        """Process one attempt at a time; only the latest can publish."""

        while True:
            with target.cond:
                while target.pending is None and not target.closed:
                    target.cond.wait()
                if target.closed:
                    return
                attempt = target.pending
                assert attempt is not None
                target.pending = None
                target.current = attempt

            # A queued attempt may already have been cancelled or superseded
            # before it ever got a CPU slot.
            if attempt.cancel_event.is_set():
                self._finish_without_publishing(target, attempt)
                attempt.terminated.set()
                continue

            attempt.status = RUNNING
            try:
                _, r1, r2 = self.resolve_target(target.job_id, target.sample_id)
                result = run_trimming(
                    source_r1=r1,
                    source_r2=r2,
                    work_dir=str(attempt.work_dir),
                    params=attempt.params,
                    is_cancelled=attempt.cancel_event.is_set,
                )
            except Exception as exc:  # noqa: BLE001 - every failure is a state
                with target.cond:
                    self._handle_failure(target, attempt, exc)
            else:
                with target.cond:
                    self._handle_success(target, attempt, result)
            finally:
                attempt.terminated.set()

    def _handle_failure(
        self, target: _TargetState, attempt: TrimAttempt, exc: Exception
    ) -> None:
        target.current = None
        attempt.finished_at = _utc_now()
        still_latest = attempt.generation == target.latest_generation
        if attempt.user_cancelled:
            attempt.status = CANCELLED
        elif not still_latest:
            attempt.status = SUPERSEDED
        else:
            attempt.status = FAILED
            attempt.error = (
                str(exc) if isinstance(exc, DemuxError) else f"unexpected error: {exc}"
            )
        # Failed / cancelled / superseded attempts never publish and leave
        # no temporary files behind.  The previously published package is
        # deliberately untouched.
        self._remove_tree(attempt.work_dir)

    def _handle_success(
        self, target: _TargetState, attempt: TrimAttempt, result: Any
    ) -> None:
        target.current = None
        attempt.finished_at = _utc_now()
        attempt.stats = dict(result.stats)
        attempt.pairs_in = result.pairs_in
        attempt.pairs_out = result.pairs_out
        attempt.zip_path = Path(result.zip_path)

        if attempt.user_cancelled or attempt.generation != target.latest_generation:
            # Lost the race while cutadapt was finishing: never publish.
            attempt.status = CANCELLED if attempt.user_cancelled else SUPERSEDED
            self._remove_tree(attempt.work_dir)
            attempt.zip_path = None
            return

        # Atomic pointer swap: the new pair + ZIP are fully complete before
        # the old package is replaced, so failures keep the old package.
        previous = target.published
        target.published = attempt
        attempt.status = COMPLETED
        if previous is not None and previous is not attempt:
            self._remove_tree(previous.work_dir)
            previous.zip_path = None

    def _finish_without_publishing(
        self, target: _TargetState, attempt: TrimAttempt
    ) -> None:
        with target.cond:
            target.current = None
            attempt.finished_at = _utc_now()
            attempt.status = (
                CANCELLED if attempt.user_cancelled else SUPERSEDED
            )
            self._remove_tree(attempt.work_dir)

    # ------------------------------------------------------------- cancel

    def cancel(self, attempt_id: str) -> TrimAttempt | None:
        """Request cancellation and wait for process exit + cleanup."""

        with self._lock:
            attempt = self._attempts.get(attempt_id)
        if attempt is None:
            return None
        if attempt.status in TERMINAL_STATES:
            return attempt
        attempt.user_cancelled = True
        attempt.cancel_event.set()
        if not attempt.terminated.wait(timeout=CANCEL_WAIT_TIMEOUT):
            raise TimeoutError(
                f"trimming attempt {attempt_id} did not stop within "
                f"{CANCEL_WAIT_TIMEOUT}s"
            )
        return attempt

    # ------------------------------------------------------------- queries

    def get_attempt(self, attempt_id: str) -> tuple[TrimAttempt, bool] | None:
        with self._lock:
            attempt = self._attempts.get(attempt_id)
            if attempt is None:
                return None
            target = self._targets.get(
                self.target_key(attempt.job_id, attempt.sample_id)
            )
            published = target is not None and target.published is attempt
            active = target.published if target is not None else None
        return attempt, published

    def _decorate(self, attempt: TrimAttempt) -> dict[str, Any]:
        """Snapshot plus the currently published package for the same sample."""
        target = self._targets.get(
            self.target_key(attempt.job_id, attempt.sample_id)
        )
        published = target is not None and target.published is attempt
        snap = attempt.snapshot(published=published)
        if target is not None and target.published is not None \
                and target.published is not attempt:
            snap["active_published"] = {
                "attempt_id": target.published.id,
                "params": target.published.params.as_dict(),
            }
        else:
            snap["active_published"] = None
        return snap

    def list_attempts(self) -> list[dict[str, Any]]:
        with self._lock:
            attempts = list(self._attempts.values())
        attempts.sort(key=lambda a: a.created_at, reverse=True)
        with self._lock:
            return [self._decorate(attempt) for attempt in attempts]

    def list_targets(self) -> list[dict[str, Any]]:
        """Non-empty samples of completed batches plus their trim state."""

        targets: list[dict[str, Any]] = []
        for job in self._demux.list_jobs_raw():
            if job.status != "completed":
                continue
            output_dir = job.work_dir / "output"
            for sample_id, count in (job.final_counts or {}).items():
                if count <= 0:
                    continue
                r1 = output_dir / f"{sample_id}_R1.fastq.gz"
                r2 = output_dir / f"{sample_id}_R2.fastq.gz"
                if not (r1.is_file() and r2.is_file() and r1.stat().st_size > 0
                        and r2.stat().st_size > 0):
                    continue
                key = self.target_key(job.id, sample_id)
                with self._lock:
                    target = self._targets.get(key)
                targets.append(
                    {
                        "target_id": key,
                        "job_id": job.id,
                        "sample_id": sample_id,
                        "read_count": count,
                        "state": target.snapshot() if target is not None else None,
                    }
                )
        targets.sort(key=lambda item: (item["job_id"], item["sample_id"]))
        return targets

    def open_published_zip(self, attempt_id: str):
        """Open the published ZIP for an attempt, keeping the fd alive.

        The file handle stays open for the whole response so a concurrent
        package replacement (which deletes the old directory) cannot cut an
        in-flight download short: the caller receives one complete version.
        """

        found = self.get_attempt(attempt_id)
        if found is None:
            raise KeyError("attempt not found (state is memory-only)")
        attempt, published = found
        if not published or attempt.zip_path is None or not attempt.zip_path.exists():
            raise FileNotFoundError("this attempt has no published ZIP")
        fh = open(attempt.zip_path, "rb")  # noqa: SIM115 - closed by response
        filename = (
            f"trimmed_{attempt.sample_id}_{attempt.job_id[:8]}_{attempt.id[:8]}.zip"
        )
        return fh, filename

    # ------------------------------------------------------------- helpers

    @staticmethod
    def _remove_tree(path: Path) -> None:
        shutil.rmtree(path, ignore_errors=True)
