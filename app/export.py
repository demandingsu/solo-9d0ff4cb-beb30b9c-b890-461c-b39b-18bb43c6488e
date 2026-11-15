"""Streaming, de-identified export of succeeded QC results.

An export reads the two published (already QC-filtered) gzip files of a
successful job pair by pair and writes a new pair of gzip files in which every
read pair gets a fresh identifier.  The original header comments are dropped,
the ``+`` line is reduced to a bare ``+``, and sequence/quality bytes are left
untouched.  The same new identifier is used for both ends (with ``/1`` and
``/2`` suffixes).

Every pair is also recorded in a mapping table that stores the new identifier
and both *original, complete* header lines.  The mapping is written only below
``DATA_DIR/exports/<export_id>/mapping.jsonl``; no API route serves it.

Failure and stop contract
-------------------------

* Source corruption, a mapping-table write failure, an output write failure
  or a publish failure removes the whole export directory (outputs, mapping,
  partial files) before ``failed`` becomes visible.
* A stop accepted while the export is pending/running moves it to
  ``stopping``; it can never subsequently become ``succeeded``.  ``stopped``
  is published only after the worker has exited and the directory is gone.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import shutil
import sys
import threading
import uuid
import zlib
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .jobs import (
    FAILED,
    PENDING,
    RUNNING,
    STOPPABLE_STATUSES,
    STOPPED,
    STOPPING,
    SUCCEEDED,
)
from .qc import (
    JobStopped,
    QCError,
    _GzipLineReader,
    _check_record,
    _strip_eol,
    _with_eol,
)

logger = logging.getLogger(__name__)

R1_OUTPUT_NAME = "R1.export.fastq.gz"
R2_OUTPUT_NAME = "R2.export.fastq.gz"
R1_PART_NAME = R1_OUTPUT_NAME + ".part"
R2_PART_NAME = R2_OUTPUT_NAME + ".part"
MAPPING_NAME = "mapping.jsonl"
MAPPING_PART_NAME = MAPPING_NAME + ".part"

PROGRESS_EVERY_PAIRS = 256


class ExportNotFound(Exception):
    pass


class ExportConflict(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


class ExportError(Exception):
    """A deterministic export failure with a machine-readable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


@dataclass(frozen=True)
class ExportResult:
    pairs: int


ProgressCallback = Callable[[int], None]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def new_read_id(export_id: str, pair_number: int) -> str:
    """Return the replacement identifier for one read pair.

    The export-specific prefix keeps identifiers unambiguous between exports
    and makes it obvious that a name was generated rather than copied from the
    source.  Nothing from the original read name is used.
    """
    return f"DX{export_id[:12]}_{pair_number:012d}"


class _ExportGzipWriter:
    """Gzip output that treats every create/write/finalize error as fatal."""

    def __init__(self, path: Path, side: str) -> None:
        self.side = side
        try:
            self._raw = open(path, "wb")
        except OSError as exc:
            raise ExportError(
                "write_failed", f"{side}: cannot create output file: {exc}"
            ) from exc
        try:
            # filename="" keeps the source/part file name out of the gzip
            # header; only the generated read identifiers appear in the data.
            self._gz = gzip.GzipFile(
                filename="",
                fileobj=self._raw,
                mode="wb",
                compresslevel=6,
                mtime=0,
            )
        except Exception as exc:
            self._raw.close()
            raise ExportError(
                "write_failed", f"{side}: cannot open output stream: {exc}"
            ) from exc
        self._closed = False

    def write(self, data: bytes) -> None:
        try:
            self._gz.write(data)
        except (OSError, ValueError, zlib.error) as exc:
            raise ExportError(
                "write_failed", f"{self.side}: cannot write output: {exc}"
            ) from exc

    def close(self) -> None:
        if self._closed:
            return
        try:
            self._gz.close()
            self._raw.flush()
            os.fsync(self._raw.fileno())
        except Exception as exc:
            raise ExportError(
                "write_failed",
                f"{self.side}: cannot finalize output file: {exc}",
            ) from exc
        finally:
            self._closed = True
            try:
                self._raw.close()
            except Exception:
                pass


class _MappingWriter:
    """JSON-lines mapping writer; only ever stored server-side."""

    def __init__(self, path: Path) -> None:
        try:
            self._raw = open(path, "wb")
        except OSError as exc:
            raise ExportError(
                "mapping_write_failed",
                f"cannot create mapping table: {exc}",
            ) from exc
        self._closed = False

    def write(
        self,
        pair_number: int,
        new_id: str,
        original_header_1: bytes,
        original_header_2: bytes,
    ) -> None:
        record = {
            "pair_number": pair_number,
            "new_id": new_id,
            "r1_header": _strip_eol(original_header_1).decode(
                "utf-8", "surrogateescape"
            ),
            "r2_header": _strip_eol(original_header_2).decode(
                "utf-8", "surrogateescape"
            ),
        }
        try:
            line = (
                json.dumps(
                    record, ensure_ascii=False, separators=(",", ":")
                )
                + "\n"
            ).encode("utf-8", "surrogateescape")
            self._raw.write(line)
        except (OSError, ValueError, UnicodeError) as exc:
            raise ExportError(
                "mapping_write_failed",
                f"cannot write mapping table: {exc}",
            ) from exc

    def close(self) -> None:
        if self._closed:
            return
        try:
            self._raw.flush()
            os.fsync(self._raw.fileno())
        except Exception as exc:
            raise ExportError(
                "mapping_write_failed",
                f"cannot finalize mapping table: {exc}",
            ) from exc
        finally:
            self._closed = True
            try:
                self._raw.close()
            except Exception:
                pass


def _record_bytes(
    new_id: str, suffix: bytes, sequence_line: bytes, quality_line: bytes
) -> bytes:
    """Build one de-identified FASTQ record.

    Only the header and plus line are replaced.  Sequence and quality content
    are copied byte for byte (line endings are normalized to ``\\n``).
    """
    header = b"@" + new_id.encode("ascii") + suffix + b"\n"
    return header + _with_eol(sequence_line) + b"+\n" + _with_eol(
        quality_line
    )


def run_export(
    source_r1: Path,
    source_r2: Path,
    output_r1_part: Path,
    output_r2_part: Path,
    mapping_part: Path,
    export_id: str,
    stop_event: threading.Event,
    *,
    progress_cb: ProgressCallback | None = None,
) -> ExportResult:
    """Stream a paired FASTQ export into temporary ``*.part`` files.

    The caller publishes the files only after this function returns, so a
    partially written export is never downloadable.
    """
    reader1 = reader2 = writer1 = writer2 = mapping = None
    try:
        reader1 = _GzipLineReader(source_r1, "R1", stop_event)
        reader2 = _GzipLineReader(source_r2, "R2", stop_event)
        writer1 = _ExportGzipWriter(output_r1_part, "R1")
        writer2 = _ExportGzipWriter(output_r2_part, "R2")
        mapping = _MappingWriter(mapping_part)

        pairs = 0
        while True:
            if stop_event.is_set():
                raise JobStopped()

            header1 = reader1.readline()
            header2 = reader2.readline()
            if not header1 and not header2:
                break
            if not header1 or not header2:
                missing = "R2" if header1 else "R1"
                other = "R1" if missing == "R2" else "R2"
                raise QCError(
                    "pair_count_mismatch",
                    f"{other} has more records than {missing}; "
                    "the two ends contain different numbers of reads",
                )

            pair_number = pairs + 1
            lines1 = [
                header1,
                reader1.readline(),
                reader1.readline(),
                reader1.readline(),
            ]
            lines2 = [
                header2,
                reader2.readline(),
                reader2.readline(),
                reader2.readline(),
            ]
            for side, lines in (("R1", lines1), ("R2", lines2)):
                if any(line == b"" for line in lines[1:]):
                    raise QCError(
                        "malformed_record",
                        f"{side}: record {pair_number} has fewer than four "
                        "lines (truncated FASTQ record)",
                    )

            id1, _ = _check_record(lines1, "R1", pair_number, b"/1")
            id2, _ = _check_record(lines2, "R2", pair_number, b"/2")
            if id1 != id2:
                raise QCError(
                    "pair_name_mismatch",
                    f"record {pair_number}: R1 read name "
                    f"{id1.decode('utf-8', 'replace')!r} does not match R2 "
                    f"read name {id2.decode('utf-8', 'replace')!r}",
                )

            if stop_event.is_set():
                raise JobStopped()

            new_id = new_read_id(export_id, pair_number)
            writer1.write(
                _record_bytes(new_id, b"/1", lines1[1], lines1[3])
            )
            writer2.write(
                _record_bytes(new_id, b"/2", lines2[1], lines2[3])
            )
            mapping.write(pair_number, new_id, lines1[0], lines2[0])
            pairs = pair_number

            if progress_cb is not None and (
                pair_number % PROGRESS_EVERY_PAIRS == 0
            ):
                progress_cb(pair_number)

        if stop_event.is_set():
            raise JobStopped()

        # Finalize explicitly so a close/fsync failure is reported as a failed
        # export rather than silently publishing an incomplete file.
        writer1.close()
        writer1 = None
        writer2.close()
        writer2 = None
        mapping.close()
        mapping = None
        return ExportResult(pairs=pairs)
    finally:
        # If the body already raised, keep that error: a secondary close
        # failure must not mask the real reason.
        failed = sys.exc_info()[0] is not None
        for reader in (reader1, reader2):
            if reader is not None:
                reader.close()
        for closer in (mapping, writer1, writer2):
            if closer is not None:
                try:
                    closer.close()
                except Exception:
                    if not failed:
                        raise


def _publish_export(export: "ExportJob") -> None:
    """Rename the completed temporary files into their download names.

    External callers never observe these names until the status is set to
    ``succeeded`` under the manager lock, and the caller removes the entire
    directory if any rename fails.
    """
    os.replace(export.mapping_part, export.mapping_path)
    os.replace(export.part_r1, export.output_r1)
    os.replace(export.part_r2, export.output_r2)


@dataclass
class ExportJob:
    export_id: str
    source_job_id: str
    source_r1: Path
    source_r2: Path
    directory: Path
    status: str = PENDING
    created_at: str = field(default_factory=_now)
    started_at: str | None = None
    finished_at: str | None = None
    pairs_processed: int = 0
    pairs_exported: int | None = None
    error_code: str | None = None
    error_message: str | None = None
    stop_event: threading.Event = field(default_factory=threading.Event)
    thread: threading.Thread | None = None

    # -- paths -----------------------------------------------------------
    @property
    def part_r1(self) -> Path:
        return self.directory / R1_PART_NAME

    @property
    def part_r2(self) -> Path:
        return self.directory / R2_PART_NAME

    @property
    def mapping_part(self) -> Path:
        return self.directory / MAPPING_PART_NAME

    @property
    def output_r1(self) -> Path:
        return self.directory / R1_OUTPUT_NAME

    @property
    def output_r2(self) -> Path:
        return self.directory / R2_OUTPUT_NAME

    @property
    def mapping_path(self) -> Path:
        return self.directory / MAPPING_NAME

    @property
    def downloads_ready(self) -> bool:
        return (
            self.status == SUCCEEDED
            and self.output_r1.is_file()
            and self.output_r2.is_file()
        )

    def to_dict(self) -> dict[str, Any]:
        total = self.pairs_exported
        result: dict[str, Any] = {
            "export_id": self.export_id,
            "job_id": self.source_job_id,
            "source_job_id": self.source_job_id,
            "status": self.status,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "pairs_processed": self.pairs_processed,
            # The exported pair count; null until the export has succeeded.
            "pairs_exported": total,
            "total_pairs": total,
            "pair_count": total,
            "error": (
                {"code": self.error_code, "message": self.error_message}
                if self.error_code
                else None
            ),
            "downloads": None,
        }
        if self.downloads_ready:
            result["downloads"] = {
                "r1": f"/exports/{self.export_id}/download/r1",
                "r2": f"/exports/{self.export_id}/download/r2",
            }
        return result


class ExportManager:
    """In-process manager for de-identified exports."""

    def __init__(self, data_dir: Path) -> None:
        self.data_dir = Path(data_dir)
        self.exports_dir = self.data_dir / "exports"
        self.exports_dir.mkdir(parents=True, exist_ok=True)
        self._exports: dict[str, ExportJob] = {}
        self._lock = threading.Lock()
        self._remove_stale_directories()

    # -- lifecycle -------------------------------------------------------
    def _remove_stale_directories(self) -> None:
        """Remove leftovers of exports that never finished.

        Export state lives in memory only, so after a restart no export can be
        queried or downloaded.  Completed directories are nevertheless kept:
        the mapping table is the long-lived, internal traceability record and
        must survive restarts.  Incomplete directories (crash, killed process)
        are discarded because they contain no usable mapping.
        """
        for child in self.exports_dir.iterdir():
            try:
                if child.is_dir():
                    complete = (
                        (child / R1_OUTPUT_NAME).is_file()
                        and (child / R2_OUTPUT_NAME).is_file()
                        and (child / MAPPING_NAME).is_file()
                    )
                    if not complete:
                        shutil.rmtree(child, ignore_errors=True)
                else:
                    child.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort
                logger.warning("could not remove stale export path %s", child)

    def create_export(
        self, source_job_id: str, source_r1: Path, source_r2: Path
    ) -> ExportJob:
        export_id = uuid.uuid4().hex
        directory = self.exports_dir / export_id
        directory.mkdir(parents=True)
        export = ExportJob(
            export_id=export_id,
            source_job_id=source_job_id,
            source_r1=Path(source_r1),
            source_r2=Path(source_r2),
            directory=directory,
        )
        with self._lock:
            self._exports[export_id] = export
        return export

    def abort(self, export_id: str) -> None:
        with self._lock:
            export = self._exports.pop(export_id, None)
        if export is not None:
            shutil.rmtree(export.directory, ignore_errors=True)

    def start(self, export_id: str) -> None:
        export = self.get(export_id)
        thread = threading.Thread(
            target=self._run_export,
            args=(export,),
            name=f"export-{export_id[:8]}",
            daemon=True,
        )
        with self._lock:
            if export.status != PENDING:
                raise ExportConflict("export has already been started")
            export.thread = thread
        thread.start()

    def shutdown(self) -> None:
        with self._lock:
            exports = list(self._exports.values())
        for export in exports:
            export.stop_event.set()
        for export in exports:
            thread = export.thread
            if thread is not None and thread.is_alive():
                thread.join(timeout=5)

    # -- queries ---------------------------------------------------------
    def get(self, export_id: str) -> ExportJob:
        with self._lock:
            export = self._exports.get(export_id)
        if export is None:
            raise ExportNotFound(export_id)
        return export

    def list_exports(
        self, source_job_id: str | None = None
    ) -> list[ExportJob]:
        with self._lock:
            exports = list(self._exports.values())
        if source_job_id is not None:
            exports = [
                export
                for export in exports
                if export.source_job_id == source_job_id
            ]
        return sorted(exports, key=lambda e: e.created_at, reverse=True)

    def snapshot(self, export_id: str) -> dict[str, Any]:
        with self._lock:
            export = self._exports.get(export_id)
            if export is None:
                raise ExportNotFound(export_id)
            return export.to_dict()

    def request_stop(self, export_id: str) -> dict[str, Any]:
        with self._lock:
            export = self._exports.get(export_id)
            if export is None:
                raise ExportNotFound(export_id)
            if export.status in STOPPABLE_STATUSES:
                export.status = STOPPING
                export.stop_event.set()
            elif export.status == STOPPING:
                export.stop_event.set()
            elif export.status == STOPPED:
                pass  # idempotent
            else:  # succeeded / failed
                raise ExportConflict(
                    f"export is already {export.status}; it can no longer "
                    "be stopped"
                )
            return export.to_dict()

    def output_paths(self, export_id: str) -> tuple[Path, Path]:
        with self._lock:
            export = self._exports.get(export_id)
            if export is None:
                raise ExportNotFound(export_id)
            if export.status != SUCCEEDED:
                raise ExportConflict(
                    f"export status is {export.status!r}; downloads are "
                    "only available after successful completion"
                )
            if not (
                export.output_r1.is_file() and export.output_r2.is_file()
            ):
                raise ExportConflict(
                    "export files are no longer available"
                )
            return export.output_r1, export.output_r2

    # -- worker ----------------------------------------------------------
    def _run_export(self, export: ExportJob) -> None:
        with self._lock:
            if export.status == STOPPING or export.stop_event.is_set():
                stop_early = True
            elif export.status == PENDING:
                export.status = RUNNING
                export.started_at = _now()
                stop_early = False
            else:  # pragma: no cover - defensive
                return
        if stop_early:
            self._finish_stopped(export)
            return

        try:
            result = run_export(
                export.source_r1,
                export.source_r2,
                export.part_r1,
                export.part_r2,
                export.mapping_part,
                export.export_id,
                export.stop_event,
                progress_cb=lambda pairs: self._set_progress(export, pairs),
            )
        except (JobStopped, QCError, ExportError) as exc:
            self._finish_failed_or_stopped(export, exc)
        except Exception as exc:  # pragma: no cover - unexpected bug
            logger.exception("unexpected error in export %s", export.export_id)
            self._finish_failed_or_stopped(
                export, ExportError("internal_error", str(exc))
            )
        else:
            self._finish_completed(export, result)

    def _set_progress(self, export: ExportJob, pairs: int) -> None:
        with self._lock:
            if pairs > export.pairs_processed:
                export.pairs_processed = pairs

    def _finish_completed(
        self, export: ExportJob, result: ExportResult
    ) -> None:
        stop_wins = False
        publish_error: Exception | None = None
        with self._lock:
            if export.status == STOPPING or export.stop_event.is_set():
                stop_wins = True
            else:
                try:
                    _publish_export(export)
                except Exception as exc:  # noqa: BLE001 - reported as failure
                    publish_error = exc
                else:
                    export.pairs_processed = result.pairs
                    export.pairs_exported = result.pairs
                    export.finished_at = _now()
                    export.status = SUCCEEDED
        if stop_wins:
            self._finish_stopped(export)
            return
        if publish_error is not None:
            logger.error(
                "could not publish export %s: %s",
                export.export_id,
                publish_error,
            )
            self._finish_failed_or_stopped(
                export,
                ExportError("publish_failed", str(publish_error)),
            )

    def _finish_failed_or_stopped(
        self, export: ExportJob, exc: BaseException
    ) -> None:
        # Clean up *before* exposing the terminal state, so a query that sees
        # "failed" or "stopped" is guaranteed that no artifacts remain.
        self._cleanup_all(export)
        with self._lock:
            export.pairs_exported = None
            export.finished_at = _now()
            if isinstance(exc, JobStopped) or (
                export.status == STOPPING or export.stop_event.is_set()
            ):
                # A stop accepted at any point wins over failure.
                export.error_code = None
                export.error_message = None
                export.status = STOPPED
            elif isinstance(exc, (QCError, ExportError)):
                export.error_code = exc.code
                export.error_message = exc.message
                export.status = FAILED
            else:  # pragma: no cover - defensive
                export.error_code = "internal_error"
                export.error_message = str(exc)
                export.status = FAILED

    def _finish_stopped(self, export: ExportJob) -> None:
        # Temporary inputs/partial outputs and any published export artifacts
        # are removed *before* the job is reported as stopped.
        self._cleanup_all(export)
        with self._lock:
            export.pairs_exported = None
            export.error_code = None
            export.error_message = None
            export.finished_at = _now()
            # Only publish ``stopped`` if the export has not already reached a
            # terminal state (defensive against repeated calls).
            if export.status not in (STOPPED, SUCCEEDED, FAILED):
                export.status = STOPPED

    def _cleanup_all(self, export: ExportJob) -> None:
        shutil.rmtree(export.directory, ignore_errors=True)
