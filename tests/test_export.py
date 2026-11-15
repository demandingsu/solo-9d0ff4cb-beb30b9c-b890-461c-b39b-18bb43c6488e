from __future__ import annotations

import gzip
import io
import json
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.main import create_app


def fastq_record(name: str, seq: str, qual: str) -> str:
    return f"@{name}\n{seq}\n+\n{qual}\n"


def gz_bytes(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(text.encode("ascii"))
    return buffer.getvalue()


def gz_pair(r1: str, r2: str, *, truncate: str | None = None):
    data1 = gz_bytes(r1)
    data2 = gz_bytes(r2)
    if truncate == "r1":
        data1 = data1[:-6]
    elif truncate == "r2":
        data2 = data2[:-6]
    return data1, data2


@pytest.fixture()
def client(tmp_path: Path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as test_client:
        yield test_client


def wait_for(client: TestClient, url: str, statuses, timeout=15.0) -> dict:
    deadline = time.monotonic() + timeout
    payload = {}
    while time.monotonic() < deadline:
        payload = client.get(url).json()
        if payload["status"] in statuses:
            return payload
        time.sleep(0.01)
    raise AssertionError(f"{url} did not reach {statuses}: {payload}")


HIGH = "I" * 4


def run_job(client, r1: str, r2: str, threshold: int = 20) -> dict:
    response = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", gz_bytes(r1), "application/gzip"),
            "r2": ("R2.fastq.gz", gz_bytes(r2), "application/gzip"),
        },
        data={"min_quality": str(threshold)},
    )
    assert response.status_code == 202, response.text
    job_id = response.json()["job_id"]
    done = wait_for(client, f"/jobs/{job_id}", {"succeeded", "failed"})
    assert done["status"] == "succeeded", done
    return done


def start_export(client, job_id: str) -> dict:
    response = client.post(f"/jobs/{job_id}/exports")
    assert response.status_code == 202, response.text
    return response.json()


def finish_export(client, export_id: str) -> dict:
    return wait_for(
        client, f"/exports/{export_id}", {"succeeded", "failed", "stopped"}
    )


def read_gz(data: bytes) -> str:
    return gzip.decompress(data).decode("ascii")


def parse_records(text: str) -> list[tuple[str, str, str, str]]:
    lines = text.splitlines()
    assert len(lines) % 4 == 0
    return [tuple(lines[i : i + 4]) for i in range(0, len(lines), 4)]


def submit_with_comments(client) -> tuple[dict, str, str]:
    r1 = (
        fastq_record("old1/1 first comment", "ACGT", "IIII")
        + fastq_record("old2/1", "TTTT", "IIII")
        + fastq_record("old3/1", "GGGG", "IIII")
    )
    r2 = (
        fastq_record("old1/2 second comment", "TGCA", "IIII")
        + fastq_record("old2/2", "CCCC", "IIII")
        + fastq_record("old3/2", "AAAA", "IIII")
    )
    # Threshold below every quality -> all pairs kept, order preserved.
    job = run_job(client, r1, r2, threshold=20)
    return job, r1, r2


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_export_replaces_identifiers_and_keeps_payload(client):
    job, r1_source, r2_source = submit_with_comments(client)
    created = start_export(client, job["job_id"])
    assert created["status"] in {"pending", "running", "succeeded"}
    assert created["job_id"] == job["job_id"]

    done = finish_export(client, created["export_id"])
    assert done["status"] == "succeeded", done
    assert done["pairs_exported"] == 3
    assert done["total_pairs"] == 3
    assert done["error"] is None
    assert done["downloads"] == {
        "r1": f"/exports/{created['export_id']}/download/r1",
        "r2": f"/exports/{created['export_id']}/download/r2",
    }

    file1 = client.get(done["downloads"]["r1"])
    file2 = client.get(done["downloads"]["r2"])
    assert file1.status_code == 200 and file2.status_code == 200
    assert file1.headers["content-type"] == "application/gzip"
    assert "deidentified" in file1.headers["content-disposition"]

    records1 = parse_records(read_gz(file1.content))
    records2 = parse_records(read_gz(file2.content))
    source1 = parse_records(r1_source)
    source2 = parse_records(r2_source)
    assert len(records1) == len(records2) == 3

    new_ids = []
    for index, (rec1, rec2) in enumerate(zip(records1, records2)):
        header1, seq1, plus1, qual1 = rec1
        header2, seq2, plus2, qual2 = rec2
        assert header1.startswith("@") and header2.startswith("@")
        assert header1.endswith("/1") and header2.endswith("/2")
        # Exactly the same identifier on both ends.
        assert header1[:-2] == header2[:-2]
        # No original name or comment leaks into either header.
        assert "old" not in header1 and "comment" not in header1
        assert "old" not in header2 and "comment" not in header2
        # Plus line carries no name.
        assert plus1 == "+" and plus2 == "+"
        # Sequence and quality unchanged.
        assert seq1 == source1[index][1] and qual1 == source1[index][3]
        assert seq2 == source2[index][1] and qual2 == source2[index][3]
        new_ids.append(header1[:-2])

    # Same pair -> same new id; different pairs get different ids; order kept.
    assert len(set(new_ids)) == 3
    records_again = parse_records(
        read_gz(client.get(done["downloads"]["r1"]).content)
    )
    assert [r[0] for r in records_again] == [r[0] for r in records1]


def test_export_keeps_only_pairs_kept_by_qc(client):
    r1 = (
        fastq_record("keep/1", "ACGT", HIGH)
        + fastq_record("drop/1", "ACGT", "III#")
        + fastq_record("keep2/1", "TTTT", HIGH)
    )
    r2 = (
        fastq_record("keep/2", "GGGG", HIGH)
        + fastq_record("drop/2", "GGGG", HIGH)
        + fastq_record("keep2/2", "CCCC", HIGH)
    )
    job = run_job(client, r1, r2, threshold=20)
    assert job["kept_pairs"] == 2
    export = finish_export(client, start_export(client, job["job_id"])["export_id"])
    assert export["pairs_exported"] == 2
    records = parse_records(
        read_gz(client.get(export["downloads"]["r1"]).content)
    )
    assert [r[1] for r in records] == ["ACGT", "TTTT"]


def test_export_mapping_is_server_side_only(client):
    job, r1_source, r2_source = submit_with_comments(client)
    export = finish_export(client, start_export(client, job["job_id"])["export_id"])
    assert export["status"] == "succeeded"

    data_dir = client.app.state.manager.data_dir
    mapping_path = data_dir / "exports" / export["export_id"] / "mapping.jsonl"
    assert mapping_path.is_file()
    rows = [
        json.loads(line) for line in mapping_path.read_text("utf-8").splitlines()
    ]
    assert [row["pair_number"] for row in rows] == [1, 2, 3]
    source_headers1 = [r[0] for r in parse_records(r1_source)]
    source_headers2 = [r[0] for r in parse_records(r2_source)]
    for index, row in enumerate(rows):
        assert row["r1_header"] == source_headers1[index]
        assert row["r2_header"] == source_headers2[index]
        assert row["r1_header"].startswith("@old")
        assert row["r2_header"].startswith("@old")
        # The mapping's new id matches the exported header for that pair.
        assert row["new_id"].endswith(f"{index + 1:012d}")

    records = parse_records(
        read_gz(client.get(export["downloads"]["r1"]).content)
    )
    assert [r[0][1:-2] for r in records] == [row["new_id"] for row in rows]

    # No API route exposes the mapping table.
    for url in (
        f"/exports/{export['export_id']}/mapping",
        f"/exports/{export['export_id']}/download/mapping",
        f"/exports/{export['export_id']}/download/mapping.jsonl",
    ):
        assert client.get(url).status_code == 404


def test_export_requires_successful_job(client):
    r1 = fastq_record("a/1", "ACGT", "III#")
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job = run_job(client, r1, r2, threshold=20)
    assert job["kept_pairs"] == 0

    failed = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", gz_bytes(fastq_record("a/1", "A", "I")), "x"),
            "r2": (
                "R2.fastq.gz",
                gz_bytes(fastq_record("a/2", "A", "I"))[:-6],
                "x",
            ),
        },
        data={"min_quality": "20"},
    ).json()
    wait_for(client, f"/jobs/{failed['job_id']}", {"failed"})
    assert client.post(f"/jobs/{failed['job_id']}/exports").status_code == 409
    assert client.post("/jobs/nope/exports").status_code == 404
    # Unknown export.
    assert client.get("/exports/nope").status_code == 404
    assert client.post("/exports/nope/stop").status_code == 404
    assert client.get("/exports/nope/download/r1").status_code == 404


def test_export_by_body(client):
    job, _, _ = submit_with_comments(client)
    response = client.post("/exports", json={"job_id": job["job_id"]})
    assert response.status_code == 202, response.text
    export = finish_export(client, response.json()["export_id"])
    assert export["status"] == "succeeded"


# ---------------------------------------------------------------------------
# Failure paths: everything is cleaned up
# ---------------------------------------------------------------------------


def _assert_no_artifacts(client, export_id):
    directory = client.app.state.manager.data_dir / "exports" / export_id
    assert not directory.exists()


def test_corrupt_source_fails_and_cleans(client, monkeypatch):
    import app.export as export_module

    job, _, _ = submit_with_comments(client)
    export_id = start_export(client, job["job_id"])["export_id"]
    finish_export(client, export_id)

    # Corrupt the published R2 after the successful QC job and export again.
    manager = client.app.state.manager
    source_path = manager.jobs_dir / job["job_id"] / "R2.kept.fastq.gz"
    data = source_path.read_bytes()
    source_path.write_bytes(data[:-6])

    second = start_export(client, job["job_id"])
    done = finish_export(client, second["export_id"])
    assert done["status"] == "failed", done
    assert done["error"]["code"] == "gzip_corrupt", done
    assert done["pairs_exported"] is None
    assert done["downloads"] is None
    _assert_no_artifacts(client, second["export_id"])
    assert client.get(f"/exports/{second['export_id']}/download/r1").status_code == 409
    # The earlier successful export and its mapping are untouched.
    assert not second["export_id"] == export_id


def test_mapping_write_failure_cleans_everything(monkeypatch, tmp_path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as client:
        job, _, _ = submit_with_comments(client)
        job_dir = app.state.manager.jobs_dir / job["job_id"]
        before = sorted(p.name for p in job_dir.iterdir())

        import app.export as export_module

        def broken_write(self, *args, **kwargs):
            raise export_module.ExportError(
                "mapping_write_failed", "disk full (simulated)"
            )

        monkeypatch.setattr(export_module._MappingWriter, "write", broken_write)
        export = start_export(client, job["job_id"])
        done = finish_export(client, export["export_id"])
        assert done["status"] == "failed", done
        assert done["error"]["code"] == "mapping_write_failed"
        assert not (
            app.state.manager.data_dir / "exports" / export["export_id"]
        ).exists()
        # Source job outputs are untouched by an export failure.
        assert sorted(p.name for p in job_dir.iterdir()) == before


def test_output_write_failure_cleans_everything(monkeypatch, tmp_path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as client:
        job, _, _ = submit_with_comments(client)

        import app.export as export_module

        def broken_write(self, data):
            raise export_module.ExportError(
                "write_failed", "no space left (simulated)"
            )

        monkeypatch.setattr(
            export_module._ExportGzipWriter, "write", broken_write
        )
        export = start_export(client, job["job_id"])
        done = finish_export(client, export["export_id"])
        assert done["status"] == "failed", done
        assert done["error"]["code"] == "write_failed"
        assert not (
            app.state.manager.data_dir / "exports" / export["export_id"]
        ).exists()


def test_publish_failure_cleans_everything(monkeypatch, tmp_path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as client:
        job, _, _ = submit_with_comments(client)

        import app.export as export_module

        def broken_publish(export):
            raise OSError("rename failed (simulated)")

        monkeypatch.setattr(export_module, "_publish_export", broken_publish)
        export = start_export(client, job["job_id"])
        done = finish_export(client, export["export_id"])
        assert done["status"] == "failed", done
        assert done["error"]["code"] == "publish_failed"
        assert not (
            app.state.manager.data_dir / "exports" / export["export_id"]
        ).exists()
        assert client.get(f"/exports/{export['export_id']}/download/r1").status_code == 409


def test_mapping_table_is_not_published_as_download(client):
    job, _, _ = submit_with_comments(client)
    export = finish_export(client, start_export(client, job["job_id"])["export_id"])
    directory = client.app.state.manager.data_dir / "exports" / export["export_id"]
    names = sorted(p.name for p in directory.iterdir())
    assert names == ["R1.export.fastq.gz", "R2.export.fastq.gz", "mapping.jsonl"]
    assert not any(name.endswith(".part") for name in names)


# ---------------------------------------------------------------------------
# Stop semantics
# ---------------------------------------------------------------------------


def _blocked_export(client, monkeypatch):
    import threading

    import app.export as export_module

    entered = threading.Event()
    release = threading.Event()

    def blocking_run_export(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return export_module.ExportResult(pairs=1)

    # Block before any file is written to exercise cleanup of a claimed
    # export directory.
    monkeypatch.setattr(export_module, "run_export", blocking_run_export)
    job, _, _ = submit_with_comments(client)
    export = start_export(client, job["job_id"])
    assert entered.wait(10)
    return export, release


def test_stop_waits_for_worker_and_cleans_before_stopped(client, monkeypatch):
    export, release = _blocked_export(client, monkeypatch)
    response = client.post(f"/exports/{export['export_id']}/stop")
    assert response.status_code == 200
    assert response.json()["status"] == "stopping"

    for _ in range(20):
        assert (
            client.get(f"/exports/{export['export_id']}").json()["status"]
            == "stopping"
        )
        time.sleep(0.005)

    release.set()
    done = finish_export(client, export["export_id"])
    assert done["status"] == "stopped", done
    assert done["pairs_exported"] is None
    directory = (
        client.app.state.manager.data_dir / "exports" / export["export_id"]
    )
    assert not directory.exists()
    assert client.get(
        f"/exports/{export['export_id']}/download/r1"
    ).status_code == 409

    # Idempotent after stopped.
    again = client.post(f"/exports/{export['export_id']}/stop")
    assert again.status_code == 200
    assert again.json()["status"] == "stopped"


def test_stop_after_success_rejected(client):
    job, _, _ = submit_with_comments(client)
    export = finish_export(client, start_export(client, job["job_id"])["export_id"])
    assert export["status"] == "succeeded"
    assert (
        client.post(f"/exports/{export['export_id']}/stop").status_code == 409
    )
    assert (
        client.get(f"/exports/{export['export_id']}/download/r1").status_code
        == 200
    )


def test_stop_wins_over_failure(client, monkeypatch):
    import threading

    import app.export as export_module

    entered = threading.Event()
    release = threading.Event()

    def failing_after_stop(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        stop_event = kwargs.get("stop_event") or args[6]
        if stop_event.is_set():
            from app.qc import JobStopped

            raise JobStopped()
        raise export_module.ExportError("write_failed", "simulated")

    monkeypatch.setattr(export_module, "run_export", failing_after_stop)
    job, _, _ = submit_with_comments(client)
    export = start_export(client, job["job_id"])
    assert entered.wait(10)
    assert client.post(f"/exports/{export['export_id']}/stop").status_code == 200
    release.set()
    done = finish_export(client, export["export_id"])
    assert done["status"] == "stopped", done
    assert done["error"] is None


def test_jobs_stay_responsive_while_export_runs(client, monkeypatch):
    import threading

    import app.export as export_module

    entered = threading.Event()
    release = threading.Event()

    def blocking_run_export(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return export_module.ExportResult(pairs=1)

    monkeypatch.setattr(export_module, "run_export", blocking_run_export)
    job, _, _ = submit_with_comments(client)
    export = start_export(client, job["job_id"])
    assert entered.wait(10)

    fresh = run_job(
        client,
        fastq_record("new/1", "ACGT", HIGH),
        fastq_record("new/2", "ACGT", HIGH),
    )
    assert fresh["status"] == "succeeded"
    assert client.get("/jobs").status_code == 200
    assert client.get(f"/jobs/{job['job_id']}").status_code == 200

    stop = client.post(f"/jobs/{job['job_id']}/stop")
    assert stop.status_code == 409  # source job succeeded
    assert client.get(f"/exports/{export['export_id']}").json()["status"] == "running"

    assert client.post(f"/exports/{export['export_id']}/stop").status_code == 200
    release.set()
    assert finish_export(client, export["export_id"])["status"] == "stopped"
    assert client.get("/health").json() == {"status": "ok"}


# ---------------------------------------------------------------------------
# Memory bounds
# ---------------------------------------------------------------------------


def test_export_streams_instead_of_loading_whole_file(monkeypatch, tmp_path):
    from app import export as export_module
    from app import qc

    calls = {"readline": 0}
    original = qc._GzipLineReader.readline

    def counting(self):
        calls["readline"] += 1
        return original(self)

    monkeypatch.setattr(qc._GzipLineReader, "readline", counting)

    total = 500
    input1 = tmp_path / "R1.fastq.gz"
    input2 = tmp_path / "R2.fastq.gz"
    input1.write_bytes(
        gz_bytes("".join(fastq_record(f"r{i}/1", "A", "I") for i in range(total)))
    )
    input2.write_bytes(
        gz_bytes("".join(fastq_record(f"r{i}/2", "T", "I") for i in range(total)))
    )
    import threading

    result = export_module.run_export(
        input1,
        input2,
        tmp_path / "o1.part",
        tmp_path / "o2.part",
        tmp_path / "mapping.part",
        "abcdef012345",
        threading.Event(),
    )
    assert result.pairs == total
    assert calls["readline"] == total * 4 * 2 + 2
    assert (tmp_path / "mapping.part").read_text().count("\n") == total


def test_large_export_succeeds_and_reports_pairs(client):
    n = 3000
    r1 = "".join(fastq_record(f"r{i}/1", "ACGT", HIGH) for i in range(n))
    r2 = "".join(fastq_record(f"r{i}/2", "TGCA", HIGH) for i in range(n))
    job = run_job(client, r1, r2)
    export = finish_export(client, start_export(client, job["job_id"])["export_id"])
    assert export["status"] == "succeeded"
    assert export["pairs_exported"] == n
    mapping = (
        client.app.state.manager.data_dir
        / "exports"
        / export["export_id"]
        / "mapping.jsonl"
    )
    assert sum(1 for _ in mapping.open()) == n


def test_docs_expose_export_endpoints(client):
    schema = client.get("/openapi.json").json()
    paths = schema["paths"]
    assert "/jobs/{job_id}/exports" in paths
    assert "/exports/{export_id}" in paths
    assert "/exports/{export_id}/stop" in paths
    assert "/exports/{export_id}/download/{side}" in paths
    assert client.get("/docs").status_code == 200
