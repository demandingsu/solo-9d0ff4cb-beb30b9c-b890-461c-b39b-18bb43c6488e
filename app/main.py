"""FastAPI application for paired-end FASTQ quality control."""

from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, RedirectResponse
from pydantic import BaseModel, Field

from . import __version__
from .export import (
    ExportConflict,
    ExportManager,
    ExportNotFound,
)
from .jobs import (
    JobConflict,
    JobManager,
    JobNotFound,
)
from .qc import MAX_PHRED, MIN_PHRED

logger = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path(os.environ.get("QC_DATA_DIR", "./data"))
UPLOAD_CHUNK_BYTES = 1 << 20
MAX_UPLOAD_BYTES = int(
    os.environ.get("QC_MAX_UPLOAD_BYTES", str(20 * 1024**3))
)


def _manager(request: Request) -> JobManager:
    return request.app.state.manager


def _exports(request: Request) -> ExportManager:
    return request.app.state.exports


class ExportRequest(BaseModel):
    job_id: str = Field(
        ...,
        min_length=1,
        description="Job id of a successfully completed QC job",
    )


async def _save_upload(
    upload: UploadFile, destination: Path, label: str
) -> int:
    """Stream one uploaded file to disk without buffering it in memory."""
    total = 0
    try:
        with open(destination, "wb") as handle:
            while True:
                chunk = await upload.read(UPLOAD_CHUNK_BYTES)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_UPLOAD_BYTES:
                    raise HTTPException(
                        status_code=413,
                        detail=(
                            f"{label} exceeds the maximum upload size of "
                            f"{MAX_UPLOAD_BYTES} bytes"
                        ),
                    )
                handle.write(chunk)
    finally:
        await upload.close()
    return total


def _parse_threshold(
    min_quality: str | None,
    quality_threshold: str | None,
    threshold: str | None,
) -> int:
    raw = None
    for candidate in (min_quality, quality_threshold, threshold):
        if candidate is not None and str(candidate).strip() != "":
            raw = str(candidate).strip()
            break
    if raw is None:
        raise HTTPException(
            status_code=422,
            detail=(
                "missing quality threshold: send 'min_quality' as a form "
                f"field (integer {MIN_PHRED}..{MAX_PHRED}, Phred+33)"
            ),
        )
    try:
        value = int(raw, 10)
    except ValueError:
        raise HTTPException(
            status_code=422,
            detail=f"quality threshold {raw!r} is not an integer",
        ) from None
    if not (MIN_PHRED <= value <= MAX_PHRED):
        raise HTTPException(
            status_code=422,
            detail=(
                f"quality threshold must be between {MIN_PHRED} and "
                f"{MAX_PHRED} (Phred+33)"
            ),
        )
    return value


def create_app(data_dir: Path | str | None = None) -> FastAPI:
    resolved_data_dir = Path(
        data_dir
        if data_dir is not None
        else os.environ.get("QC_DATA_DIR", str(DEFAULT_DATA_DIR))
    )

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        manager = JobManager(resolved_data_dir)
        exports = ExportManager(resolved_data_dir)
        app.state.manager = manager
        app.state.exports = exports
        try:
            yield
        finally:
            exports.shutdown()
            manager.shutdown()

    app = FastAPI(
        title="Paired-end FASTQ QC API",
        version=__version__,
        description=(
            "Streaming quality-control filtering for paired-end FASTQ files.\n\n"
            "Upload `R1`/`R2` `.fastq.gz` files together with a Phred+33 "
            "quality threshold. A pair is kept only when every base at both "
            "ends has quality at least the threshold; otherwise the whole "
            "pair is discarded. Pair order is preserved.\n\n"
            "A job can be polled, stopped, and (after successful completion "
            "only) downloaded.\n\n"
            "For a successfully completed QC job, a de-identified export can "
            "be submitted: it streams the published results again, replaces "
            "every read name with a new generated identifier (same identifier "
            "for both ends, with `/1` and `/2` suffixes), drops original "
            "header comments and `+` line annotations, and keeps sequence "
            "and quality unchanged. The server keeps an internal read-level "
            "mapping table under its data directory; the export API never "
            "serves it. Exports can be polled, stopped, and (after success "
            "only) downloaded."
        ),
        lifespan=lifespan,
    )

    @app.get("/", include_in_schema=False)
    async def index() -> RedirectResponse:
        return RedirectResponse(url="/docs")

    @app.get("/health", tags=["meta"])
    async def health() -> dict[str, str]:
        return {"status": "ok"}

    @app.post(
        "/jobs",
        status_code=202,
        tags=["jobs"],
        summary="Upload a read pair and start a QC job",
        responses={
            202: {"description": "Job accepted and queued"},
            422: {"description": "Invalid threshold or missing files"},
        },
    )
    async def create_job(
        request: Request,
        r1: UploadFile = File(..., description="Read 1 FASTQ.gz"),
        r2: UploadFile = File(..., description="Read 2 FASTQ.gz"),
        min_quality: str | None = Form(
            None,
            description=(
                f"Minimum Phred+33 quality ({MIN_PHRED}..{MAX_PHRED}). "
                "Pairs with any base below it are discarded."
            ),
        ),
        quality_threshold: str | None = Form(None, include_in_schema=False),
        threshold: str | None = Form(None, include_in_schema=False),
    ) -> dict:
        threshold_value = _parse_threshold(
            min_quality, quality_threshold, threshold
        )
        manager = _manager(request)
        job = manager.create_job(
            r1_filename=Path(r1.filename or "R1.fastq.gz").name,
            r2_filename=Path(r2.filename or "R2.fastq.gz").name,
            min_quality=threshold_value,
        )
        try:
            size1 = await _save_upload(r1, job.input_r1, "R1")
            size2 = await _save_upload(r2, job.input_r2, "R2")
            if size1 == 0 or size2 == 0:
                raise HTTPException(
                    status_code=400,
                    detail=(
                        "one of the uploaded files is empty: "
                        f"R1={size1} bytes, R2={size2} bytes"
                    ),
                )
            manager.start(job.job_id)
        except HTTPException:
            manager.abort_upload(job.job_id)
            raise
        except Exception:
            manager.abort_upload(job.job_id)
            raise
        return manager.snapshot(job.job_id)

    @app.get("/jobs", tags=["jobs"], summary="List jobs (newest first)")
    async def list_jobs(request: Request) -> dict:
        manager = _manager(request)
        return {"jobs": [job.to_dict() for job in manager.list_jobs()]}

    @app.get(
        "/jobs/{job_id}",
        tags=["jobs"],
        summary="Get job status and counts",
        responses={404: {"description": "Unknown job"}},
    )
    async def get_job(request: Request, job_id: str) -> dict:
        try:
            return _manager(request).snapshot(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None

    @app.post(
        "/jobs/{job_id}/stop",
        tags=["jobs"],
        summary="Request job cancellation",
        responses={
            404: {"description": "Unknown job"},
            409: {"description": "Job already finished"},
        },
    )
    async def stop_job(request: Request, job_id: str) -> dict:
        try:
            return _manager(request).request_stop(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None

    @app.get(
        "/jobs/{job_id}/download/{side}",
        tags=["jobs"],
        summary="Download a filtered gzip file (successful jobs only)",
        response_class=FileResponse,
        responses={
            404: {"description": "Unknown job"},
            409: {"description": "Job not successfully completed"},
        },
    )
    async def download(
        request: Request, job_id: str, side: str
    ) -> FileResponse:
        if side not in {"r1", "r2"}:
            raise HTTPException(
                status_code=404, detail="side must be 'r1' or 'r2'"
            )
        try:
            r1_path, r2_path = _manager(request).output_paths(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None
        path, filename = (
            (r1_path, "R1.kept.fastq.gz")
            if side == "r1"
            else (r2_path, "R2.kept.fastq.gz")
        )
        return FileResponse(
            path,
            media_type="application/gzip",
            filename=filename,
        )

    # -- de-identified exports ------------------------------------------
    @app.post(
        "/jobs/{job_id}/exports",
        status_code=202,
        tags=["exports"],
        summary="Start a de-identified export of a successful QC job",
        responses={
            202: {"description": "Export accepted and queued"},
            404: {"description": "Unknown job or export"},
            409: {"description": "Job not successfully completed"},
        },
    )
    async def create_export(request: Request, job_id: str) -> dict:
        try:
            source_r1, source_r2 = _manager(request).output_paths(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        except JobConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None

        exports = _exports(request)
        export = exports.create_export(job_id, source_r1, source_r2)
        try:
            exports.start(export.export_id)
        except Exception:
            exports.abort(export.export_id)
            raise
        return exports.snapshot(export.export_id)

    @app.post(
        "/exports",
        status_code=202,
        tags=["exports"],
        summary="Start a de-identified export by body-provided job id",
        responses={
            202: {"description": "Export accepted and queued"},
            404: {"description": "Unknown job"},
            409: {"description": "Job not successfully completed"},
            422: {"description": "Missing or invalid job_id"},
        },
    )
    async def create_export_from_body(
        request: Request, payload: ExportRequest
    ) -> dict:
        return await create_export(request, payload.job_id)

    @app.get(
        "/exports",
        tags=["exports"],
        summary="List all exports (newest first)",
    )
    async def list_exports(request: Request) -> dict:
        return {
            "exports": [
                export.to_dict()
                for export in _exports(request).list_exports()
            ]
        }

    @app.get(
        "/jobs/{job_id}/exports",
        tags=["exports"],
        summary="List exports started from a job (newest first)",
        responses={404: {"description": "Unknown job"}},
    )
    async def list_job_exports(request: Request, job_id: str) -> dict:
        try:
            _manager(request).get(job_id)
        except JobNotFound:
            raise HTTPException(
                status_code=404, detail=f"job {job_id!r} not found"
            ) from None
        return {
            "exports": [
                export.to_dict()
                for export in _exports(request).list_exports(job_id)
            ]
        }

    @app.get(
        "/exports/{export_id}",
        tags=["exports"],
        summary="Get export status and exported pair count",
        responses={404: {"description": "Unknown export"}},
    )
    async def get_export(request: Request, export_id: str) -> dict:
        try:
            return _exports(request).snapshot(export_id)
        except ExportNotFound:
            raise HTTPException(
                status_code=404, detail=f"export {export_id!r} not found"
            ) from None

    @app.post(
        "/exports/{export_id}/stop",
        tags=["exports"],
        summary="Request export cancellation",
        responses={
            404: {"description": "Unknown export"},
            409: {"description": "Export already finished"},
        },
    )
    async def stop_export(request: Request, export_id: str) -> dict:
        try:
            return _exports(request).request_stop(export_id)
        except ExportNotFound:
            raise HTTPException(
                status_code=404, detail=f"export {export_id!r} not found"
            ) from None
        except ExportConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None

    @app.get(
        "/exports/{export_id}/download/{side}",
        tags=["exports"],
        summary="Download a de-identified gzip file (successful exports only)",
        response_class=FileResponse,
        responses={
            404: {"description": "Unknown export"},
            409: {"description": "Export not successfully completed"},
        },
    )
    async def download_export(
        request: Request, export_id: str, side: str
    ) -> FileResponse:
        if side not in {"r1", "r2"}:
            raise HTTPException(
                status_code=404, detail="side must be 'r1' or 'r2'"
            )
        try:
            r1_path, r2_path = _exports(request).output_paths(export_id)
        except ExportNotFound:
            raise HTTPException(
                status_code=404, detail=f"export {export_id!r} not found"
            ) from None
        except ExportConflict as exc:
            raise HTTPException(status_code=409, detail=exc.message) from None
        path, filename = (
            (r1_path, "R1.deidentified.fastq.gz")
            if side == "r1"
            else (r2_path, "R2.deidentified.fastq.gz")
        )
        return FileResponse(
            path,
            media_type="application/gzip",
            filename=filename,
        )

    return app


app = create_app()
