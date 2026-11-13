"""FastAPI application: upload, monitor, cancel and download demux jobs."""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, StreamingResponse
from fastapi.templating import Jinja2Templates

from .errors import DemuxError, TrimmingError
from .jobs import (
    MAX_PARALLEL_JOBS,
    JobManager,
)
from .trim_jobs import TrimManager
from .trimming import validate_params

BASE_DIR = Path(__file__).resolve().parent
JOBS_ROOT = Path(os.environ.get("DEMUX_JOBS_ROOT", "/data/jobs"))
TRIM_ROOT = Path(os.environ.get("DEMUX_TRIM_ROOT", "/data/trimming"))
LANES = ("R1", "R2", "I1", "I2")
SHEET_FIELD = "sample_sheet"

@asynccontextmanager
async def lifespan(_app: FastAPI):
    # Task state lives in memory only; start every run with a clean disk.
    manager.wipe()
    trim_manager.wipe()
    yield


app = FastAPI(title="FASTQ Demultiplexer", version="1.0.0", lifespan=lifespan)
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))
manager = JobManager(JOBS_ROOT)
trim_manager = TrimManager(TRIM_ROOT, manager)


# ---------------------------------------------------------------- helpers


async def _stream_to_file(upload: UploadFile, dest: Path) -> None:
    """Save an upload chunk by chunk so whole files are never buffered."""
    try:
        with open(dest, "wb") as fh:
            while True:
                chunk = await upload.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"failed to store upload: {exc}")
    finally:
        await upload.close()


# ----------------------------------------------------------------- routes


@app.get("/", response_class=HTMLResponse)
async def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        request,
        "index.html",
        {"max_parallel": MAX_PARALLEL_JOBS},
    )


@app.get("/api/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/jobs", status_code=201)
async def create_job(
    R1: UploadFile = File(...),
    R2: UploadFile = File(...),
    I1: UploadFile = File(...),
    I2: UploadFile = File(...),
    sample_sheet: UploadFile = File(...),
) -> JSONResponse:
    job = manager.try_create()
    if job is None:
        raise HTTPException(
            status_code=409,
            detail=(
                f"already running the maximum of {MAX_PARALLEL_JOBS} batches; "
                "wait for one to finish or cancel it and retry"
            ),
        )

    uploads = {
        "R1": (R1, "R1.fastq.gz"),
        "R2": (R2, "R2.fastq.gz"),
        "I1": (I1, "I1.fastq.gz"),
        "I2": (I2, "I2.fastq.gz"),
        SHEET_FIELD: (sample_sheet, "samples.csv"),
    }
    upload_names: dict[str, str] = {}
    try:
        for _key, (upload, target_name) in uploads.items():
            if not upload.filename:
                raise HTTPException(status_code=400, detail=f"missing filename for {target_name}")
            upload_names[target_name] = upload.filename
            await _stream_to_file(upload, job.upload_dir / target_name)
    except HTTPException:
        manager.discard(job)
        raise
    except Exception as exc:  # noqa: BLE001
        manager.discard(job)
        raise HTTPException(status_code=400, detail=f"upload failed: {exc}") from exc

    manager.launch(job, upload_names)
    return JSONResponse(status_code=201, content=job.snapshot())


@app.get("/api/jobs")
async def list_jobs() -> dict[str, list[dict]]:
    return {"jobs": manager.list_jobs()}


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str) -> dict:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    return job.snapshot()


@app.post("/api/jobs/{job_id}/cancel")
async def cancel_job(job_id: str) -> dict:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    try:
        job = manager.cancel(job_id)
    except TimeoutError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    assert job is not None
    return job.snapshot()


@app.get("/api/jobs/{job_id}/download")
async def download_job(job_id: str) -> FileResponse:
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found (state is memory-only)")
    if job.status != "completed" or job.zip_path is None or not job.zip_path.exists():
        raise HTTPException(
            status_code=404,
            detail=(
                "no download available: only completed jobs publish a ZIP "
                f"(current status: {job.status})"
            ),
        )
    return FileResponse(
        job.zip_path,
        media_type="application/zip",
        filename=f"demux_{job_id}.zip",
    )


@app.exception_handler(DemuxError)
async def demux_error_handler(_request: Request, exc: DemuxError) -> JSONResponse:
    # Pipeline errors are recorded on the job; this is a defensive fallback.
    return JSONResponse(status_code=400, content={"detail": str(exc)})


# ------------------------------------------------------------- trimming


@app.get("/api/trimming/targets")
async def list_trim_targets() -> dict[str, list[dict]]:
    """Non-empty samples of completed batches, with their trimming state."""
    return {"targets": trim_manager.list_targets()}


@app.post("/api/trimming", status_code=201)
async def create_trimming(
    job_id: str = Form(...),
    sample_id: str = Form(...),
    adapter_r1: str = Form(...),
    adapter_r2: str = Form(...),
    min_length: int = Form(...),
) -> JSONResponse:
    # A resubmission for the same sample must always be accepted; validation
    # failures are 400, an unknown/empty target is 404.
    try:
        params = validate_params(adapter_r1, adapter_r2, min_length)
    except TrimmingError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    try:
        trim_manager.resolve_target(job_id, sample_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    attempt = trim_manager.submit(job_id, sample_id, params)
    return JSONResponse(status_code=201, content=trim_manager._decorate(attempt))


@app.get("/api/trimming")
async def list_trimming() -> dict[str, list[dict]]:
    return {"attempts": trim_manager.list_attempts()}


@app.get("/api/trimming/{attempt_id}")
async def get_trimming(attempt_id: str) -> dict:
    found = trim_manager.get_attempt(attempt_id)
    if found is None:
        raise HTTPException(
            status_code=404, detail="trimming attempt not found (state is memory-only)"
        )
    return trim_manager._decorate(found[0])


@app.post("/api/trimming/{attempt_id}/cancel")
def cancel_trimming(attempt_id: str) -> dict:
    # Defined as a sync route on purpose: cancellation waits for the real
    # cutadapt process group to exit and the temp directory to be removed.
    # FastAPI runs sync routes in its threadpool, so the wait never blocks
    # the event loop, other samples or the demultiplex service.
    try:
        attempt = trim_manager.cancel(attempt_id)
    except TimeoutError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    if attempt is None:
        raise HTTPException(
            status_code=404, detail="trimming attempt not found (state is memory-only)"
        )
    found = trim_manager.get_attempt(attempt_id)
    assert found is not None
    return trim_manager._decorate(found[0])


@app.get("/api/trimming/{attempt_id}/download")
async def download_trimming(attempt_id: str) -> StreamingResponse:
    try:
        fh, filename = trim_manager.open_published_zip(attempt_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=404,
            detail=(
                "no download available: only the published, fully completed "
                f"attempt serves a ZIP ({exc})"
            ),
        ) from exc

    def _iterfile():
        try:
            while True:
                chunk = fh.read(1024 * 1024)
                if not chunk:
                    break
                yield chunk
        finally:
            # The fd stays open for the whole download even if the attempt
            # is superseded and its directory removed meanwhile.
            fh.close()

    return StreamingResponse(
        _iterfile(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
