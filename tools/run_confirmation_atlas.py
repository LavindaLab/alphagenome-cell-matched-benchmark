#!/usr/bin/env python3
"""Resumable, outcome-blind Atlas+AVI retrieval for the frozen confirmation set."""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import math
import os
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

from run_atlas_metadata_check import (
    CONTEXTS,
    import_sdk,
    json_safe,
    metadata_records,
    normalize_scores,
    select_required_tracks,
)


SCORERS = {"DNASE": "DNASE", "H3K27ac": "CHIP_HISTONE"}
AVI_SCORER = "AVI_SCORE"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def load_manifest(path: Path, expected: int):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "variant_id", "chromosome", "position_1based", "reference",
            "alternate", "source_variant_id_hg19",
        }
        if not required.issubset(fields):
            raise ValueError(f"missing manifest fields: {sorted(required - fields)}")
        forbidden = {
            name for name in fields
            if any(token in name.lower() for token in ("emvar", "skew", "fdr", "pvalue", "outcome"))
        }
        if forbidden:
            raise ValueError(f"outcome-like manifest fields forbidden: {sorted(forbidden)}")
        records = list(reader)
    if len(records) != expected:
        raise ValueError(f"expected {expected} variants, found {len(records)}")
    if len({row["variant_id"] for row in records}) != len(records):
        raise ValueError("duplicate hg38 variant IDs")
    return records


def canonical_bytes(records) -> bytes:
    payload = "\n".join(json.dumps(row, sort_keys=True, separators=(",", ":")) for row in records)
    return (payload + ("\n" if payload else "")).encode()


def write_gzip_atomic(path: Path, data: bytes) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(data)
    temporary.replace(path)


def atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(value)
    temporary.replace(path)


class RateLimiter:
    def __init__(self, calls_per_minute: int):
        self.limit = calls_per_minute
        self.calls = collections.deque()
        self.lock = threading.Lock()

    def wait(self) -> None:
        while True:
            with self.lock:
                moment = time.monotonic()
                while self.calls and moment - self.calls[0] >= 60:
                    self.calls.popleft()
                if len(self.calls) < self.limit:
                    self.calls.append(moment)
                    return
                delay = max(0.01, 60 - (moment - self.calls[0]) + 0.01)
            time.sleep(delay)


class ManagedRun:
    def __init__(self, directory: Path | None):
        self.directory = directory
        if directory is not None:
            directory.mkdir(parents=True, exist_ok=True)
            self.write("RUN_PID", f"{os.getpid()}\n")
            self.write("RUN_STATUS", "running\n")
            self.stage("initializing")
            self.event("run_started")

    def write(self, name: str, value: str) -> None:
        if self.directory is not None:
            atomic_text(self.directory / name, value)

    def event(self, message: str, **details) -> None:
        if self.directory is None:
            return
        record = {"utc": now(), "message": message, **details}
        with (self.directory / "events.jsonl").open("a") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
        self.write("HEARTBEAT", record["utc"] + "\n")

    def stage(self, value: str) -> None:
        self.write("CURRENT_STAGE", value + "\n")
        self.write("HEARTBEAT", now() + "\n")

    def complete(self, report: str) -> None:
        self.write("FINAL_REPORT.md", report)
        self.stage("retrieval_complete_pending_independent_inspection")
        self.write("RUN_STATUS", "completed\n")
        self.write("COMPLETED", now() + "\n")
        self.event("run_completed")

    def fail(self, error: str) -> None:
        self.write("FINAL_REPORT.md", f"# AlphaGenome confirmation retrieval\n\nFAILED: {error}\n")
        self.stage("failed")
        self.write("RUN_STATUS", "failed\n")
        self.write("FAILED", now() + "\n")
        self.event("run_failed", error=error)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--expected-variants", type=int, required=True)
    parser.add_argument("--sdk-path", type=Path, required=True)
    parser.add_argument("--pilot-metadata", type=Path, required=True)
    parser.add_argument("--managed-run-dir", type=Path)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--requests-per-minute", type=int, default=900)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    managed = ManagedRun(args.managed_run_dir)
    try:
        records = load_manifest(args.manifest_tsv, args.expected_variants)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = args.output_dir / "cache_chunks"
        cache_dir.mkdir(parents=True, exist_ok=True)
        sdk_commit = subprocess.run(
            ["git", "-C", str(args.sdk_path), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        preflight = {
            "generated_utc": now(),
            "purpose": "outcome-blind held-out confirmation Atlas and AVI retrieval",
            "manifest_sha256": sha256(args.manifest_tsv),
            "pilot_metadata_sha256": sha256(args.pilot_metadata),
            "variants": len(records),
            "requested_scorers": [*SCORERS.values(), AVI_SCORER],
            "ontology_terms": CONTEXTS,
            "combined_base_requests": len(records),
            "chunk_size": args.chunk_size,
            "max_workers": args.max_workers,
            "requests_per_minute_cap": args.requests_per_minute,
            "retries": args.retries,
            "projected_minutes_before_retries": len(records) / args.requests_per_minute,
            "sdk_commit": sdk_commit,
            "mpra_workbook_read": False,
            "mpra_functional_outcomes_accessed": False,
            "scientific_rules_mutated": False,
        }
        atomic_text(args.output_dir / "atlas_preflight.json", json.dumps(preflight, indent=2, sort_keys=True) + "\n")
        if args.preflight_only:
            print(f"preflight=pass variants={len(records)} requests={len(records)}")
            return

        key = os.environ.get("ALPHA_GENOME_API_KEY")
        if not key:
            raise RuntimeError("ALPHA_GENOME_API_KEY is not set")
        atlas, genome, _, sdk_version = import_sdk(args.sdk_path)
        client = atlas.create(key, timeout=30)
        metadata = client.scorer_metadata()
        required_scorers = {*SCORERS.values(), AVI_SCORER}
        if not required_scorers.issubset(metadata):
            raise RuntimeError(f"missing scorers: {sorted(required_scorers - set(metadata))}")
        current_metadata = {}
        for modality, scorer in SCORERS.items():
            info = metadata[scorer]
            current_metadata[modality] = {
                "scorer_name": scorer,
                "is_signed": info.is_signed,
                "required_tracks": select_required_tracks(metadata_records(info.track_metadata), modality),
            }
        pilot_metadata = json.loads(args.pilot_metadata.read_text())
        current_pairs = {
            modality: sorted((row["ontology_curie"], row["track_index"]) for row in value["required_tracks"])
            for modality, value in current_metadata.items()
        }
        pilot_pairs = {
            modality: sorted((row["ontology_curie"], row["track_index"]) for row in pilot_metadata[modality]["required_tracks"])
            for modality in SCORERS
        }
        if current_pairs != pilot_pairs:
            raise RuntimeError(f"frozen scorer-track identifiers changed: current={current_pairs} pilot={pilot_pairs}")
        metadata_receipt = {
            "retrieved_utc": now(),
            "sdk_version": sdk_version,
            "sdk_commit": sdk_commit,
            "frozen_track_pairs_match_pilot": True,
            "track_pairs": current_pairs,
            "avi_is_signed": metadata[AVI_SCORER].is_signed,
            "avi_track_metadata": metadata_records(metadata[AVI_SCORER].track_metadata),
            "mpra_functional_outcomes_accessed": False,
        }
        atomic_text(args.output_dir / "atlas_scorer_metadata.json", json.dumps(metadata_receipt, indent=2, sort_keys=True) + "\n")
        metadata_sha = sha256(args.output_dir / "atlas_scorer_metadata.json")
        limiter = RateLimiter(args.requests_per_minute)
        started = time.monotonic()
        chunks_total = math.ceil(len(records) / args.chunk_size)
        managed.stage(f"retrieving_chunk_0_of_{chunks_total}")

        def retrieve(row):
            variant = genome.Variant(
                chromosome=row["chromosome"],
                position=int(row["position_1based"]),
                reference_bases=row["reference"],
                alternate_bases=row["alternate"],
            )
            error = ""
            for attempt in range(1, args.retries + 2):
                limiter.wait()
                try:
                    response = client.query_variant(
                        variant,
                        requested_scorers=[*SCORERS.values(), AVI_SCORER],
                        ontology_terms=list(CONTEXTS),
                    )
                    regulatory = normalize_scores(row, response, SCORERS)
                    avi_data = response[AVI_SCORER]
                    quantiles = avi_data.layers.get("quantiles")
                    avi = {
                        "scorer_name": AVI_SCORER,
                        "raw_score": json_safe(avi_data.X[0, 0]),
                        "quantile_score": json_safe(quantiles[0, 0]) if quantiles is not None else None,
                    }
                    pairs = {
                        (item["ontology_curie"], item["modality"])
                        for item in regulatory
                        if item["raw_score"] is not None and item["quantile_score"] is not None
                    }
                    avi_valid = all(
                        isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))
                        for value in (avi["raw_score"], avi["quantile_score"])
                    )
                    complete = len(pairs) == 10 and avi_valid
                    normalized = {"regulatory_scores": regulatory, "avi": avi}
                    return {
                        "variant_id_hg38": row["variant_id"],
                        "source_variant_id_hg19": row["source_variant_id_hg19"],
                        **normalized,
                        "status": "complete" if complete else "incomplete",
                        "attempts": attempt,
                        "error": "." if complete else "incomplete required scores",
                        "scorer_metadata_sha256": metadata_sha,
                        "response_sha256": hashlib.sha256(canonical_bytes([normalized])).hexdigest(),
                        "mpra_functional_outcomes_accessed": False,
                    }
                except Exception as exc:  # bounded and persisted below
                    error = " ".join(f"{type(exc).__name__}: {exc}".split())
                    if attempt <= args.retries:
                        time.sleep(min(2 ** (attempt - 1), 4))
            return {
                "variant_id_hg38": row["variant_id"],
                "source_variant_id_hg19": row["source_variant_id_hg19"],
                "regulatory_scores": [],
                "avi": None,
                "status": "incomplete",
                "attempts": args.retries + 1,
                "error": error,
                "scorer_metadata_sha256": metadata_sha,
                "response_sha256": None,
                "mpra_functional_outcomes_accessed": False,
            }

        chunk_receipts = []
        for chunk_index, start in enumerate(range(0, len(records), args.chunk_size), start=1):
            chunk_records = records[start:start + args.chunk_size]
            chunk_path = cache_dir / f"chunk_{chunk_index:04d}.jsonl.gz"
            receipt_path = cache_dir / f"chunk_{chunk_index:04d}.receipt.json"
            if receipt_path.exists() and chunk_path.exists():
                old = json.loads(receipt_path.read_text())
                if old.get("passed") and old.get("file_sha256") == sha256(chunk_path) and old.get("records") == len(chunk_records):
                    chunk_receipts.append(old)
                    managed.stage(f"retrieving_chunk_{chunk_index}_of_{chunks_total}_resumed")
                    continue
            managed.stage(f"retrieving_chunk_{chunk_index}_of_{chunks_total}")
            chunk_started = time.monotonic()
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                results = list(executor.map(retrieve, chunk_records))
            results.sort(key=lambda row: row["variant_id_hg38"])
            uncompressed = canonical_bytes(results)
            write_gzip_atomic(chunk_path, uncompressed)
            counts = collections.Counter(row["status"] for row in results)
            receipt = {
                "chunk": chunk_index,
                "completed_utc": now(),
                "records": len(results),
                "complete": counts["complete"],
                "incomplete": counts["incomplete"],
                "attempts": sum(int(row["attempts"]) for row in results),
                "elapsed_seconds": round(time.monotonic() - chunk_started, 3),
                "uncompressed_sha256": hashlib.sha256(uncompressed).hexdigest(),
                "file_sha256": sha256(chunk_path),
                "file_bytes": chunk_path.stat().st_size,
                "passed": len(results) == len(chunk_records),
                "mpra_functional_outcomes_accessed": False,
            }
            atomic_text(receipt_path, json.dumps(receipt, indent=2, sort_keys=True) + "\n")
            chunk_receipts.append(receipt)
            managed.event("chunk_completed", chunk=chunk_index, complete=counts["complete"], incomplete=counts["incomplete"])

        completed = sum(item["complete"] for item in chunk_receipts)
        incomplete = sum(item["incomplete"] for item in chunk_receipts)
        attempts = sum(item["attempts"] for item in chunk_receipts)
        tree_payload = "\n".join(f"{item['chunk']:04d}\t{item['file_sha256']}" for item in chunk_receipts) + "\n"
        summary = {
            **preflight,
            "completed_utc": now(),
            "elapsed_seconds_this_invocation": round(time.monotonic() - started, 3),
            "chunks": len(chunk_receipts),
            "complete": completed,
            "incomplete": incomplete,
            "complete_fraction": completed / len(records),
            "request_attempts": attempts,
            "cache_tree_sha256": hashlib.sha256(tree_payload.encode()).hexdigest(),
            "cache_compressed_bytes": sum(item["file_bytes"] for item in chunk_receipts),
            "frozen_track_pairs_match_pilot": True,
            "api_key_persisted": False,
        }
        atomic_text(args.output_dir / "retrieval_summary.json", json.dumps(summary, indent=2, sort_keys=True) + "\n")
        report = (
            "# AlphaGenome confirmation retrieval\n\n"
            f"- Complete: {completed}/{len(records)} ({summary['complete_fraction']:.4%})\n"
            f"- Incomplete: {incomplete}\n"
            f"- Request attempts: {attempts}\n"
            f"- Cache chunks: {len(chunk_receipts)}\n"
            f"- Compressed bytes: {summary['cache_compressed_bytes']}\n"
            "- MPRA confirmation outcomes accessed: no\n"
            "- Independent inspection: pending\n"
        )
        managed.complete(report)
        print(f"complete={completed}/{len(records)} incomplete={incomplete} attempts={attempts}")
    except Exception as exc:
        managed.fail(" ".join(f"{type(exc).__name__}: {exc}".split()))
        raise


if __name__ == "__main__":
    main()
