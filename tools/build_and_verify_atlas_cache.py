#!/usr/bin/env python3
"""Materialize and deterministically replay per-variant Atlas cache records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_lines(rows: list[dict[str, object]]) -> str:
    rows.sort(
        key=lambda item: (
            str(item["variant_id_hg38"]),
            str(item["ontology_curie"]),
            str(item["modality"]),
            str(item["track_index"]),
        )
    )
    return "\n".join(json.dumps(row, sort_keys=True) for row in rows) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("cache_dir", type=Path)
    parser.add_argument("cache_manifest_json", type=Path)
    parser.add_argument("replay_receipt_json", type=Path)
    args = parser.parse_args()

    preflight = json.loads((args.run_dir / "atlas_preflight.json").read_text())
    summary = json.loads((args.run_dir / "atlas_retrieval_summary.json").read_text())
    metadata_path = args.run_dir / "atlas_scorer_metadata.json"
    metadata_sha256 = sha256_bytes(metadata_path.read_bytes())
    score_path = args.run_dir / "atlas_scores.jsonl"
    score_rows = [json.loads(line) for line in score_path.read_text().splitlines()]
    by_variant: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in score_rows:
        by_variant[str(row["variant_id_hg38"])].append(row)
    with (args.run_dir / "atlas_variant_status.tsv").open(newline="") as handle:
        status = {
            row["variant_id_hg38"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    cache_entries = []
    for variant, rows in sorted(by_variant.items()):
        normalized = canonical_lines(rows)
        response_sha256 = sha256_bytes(normalized.encode())
        record = {
            "cache_schema_version": 1,
            "normalized_variant_hg38": variant,
            "source_variant_id_hg19": rows[0]["source_variant_id_hg19"],
            "requested_ontology_terms": preflight["ontology_terms"],
            "requested_scorers": preflight["requested_scorers"],
            "raw_and_quantile_track_rows": rows,
            "scorer_metadata_sha256": metadata_sha256,
            "request_time_utc": summary["completed_utc"],
            "request_time_scope": "run-completion timestamp; API exposes no per-response timestamp",
            "sdk_version": preflight["sdk_version"],
            "sdk_commit": preflight["sdk_commit"],
            "api_resource_version": None,
            "api_resource_version_note": "not exposed by the response",
            "response_sha256": response_sha256,
            "request_attempts": int(status[variant]["attempts"]),
            "retry_history": [],
            "mpra_functional_outcomes_accessed": False,
        }
        data = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
        filename = f"{hashlib.sha256(variant.encode()).hexdigest()}.json"
        destination = args.cache_dir / filename
        temporary = args.cache_dir / f".{filename}.tmp"
        temporary.write_bytes(data)
        temporary.replace(destination)
        cache_entries.append(
            {
                "variant_id_hg38": variant,
                "file": filename,
                "file_sha256": sha256_bytes(data),
                "response_sha256": response_sha256,
                "bytes": len(data),
            }
        )

    tree_payload = "\n".join(
        f"{entry['file']}\t{entry['file_sha256']}" for entry in cache_entries
    ) + "\n"
    cache_manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(cache_entries),
        "total_bytes": sum(int(entry["bytes"]) for entry in cache_entries),
        "tree_sha256": sha256_bytes(tree_payload.encode()),
        "source_scores_sha256": sha256_bytes(score_path.read_bytes()),
        "scorer_metadata_sha256": metadata_sha256,
        "record_fields": list(record),
        "mpra_functional_outcomes_accessed": False,
    }
    args.cache_manifest_json.write_text(
        json.dumps(cache_manifest, indent=2, sort_keys=True) + "\n"
    )

    replayed_rows = []
    record_checks = []
    for entry in cache_entries:
        data = (args.cache_dir / str(entry["file"])).read_bytes()
        record = json.loads(data)
        normalized = canonical_lines(record["raw_and_quantile_track_rows"])
        record_checks.append(
            sha256_bytes(data) == entry["file_sha256"]
            and sha256_bytes(normalized.encode()) == record["response_sha256"]
        )
        replayed_rows.extend(record["raw_and_quantile_track_rows"])
    replayed = canonical_lines(replayed_rows).encode()
    source = score_path.read_bytes()
    replay_receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "cache_records": len(cache_entries),
        "all_record_and_response_checksums_valid": all(record_checks),
        "source_scores_sha256": sha256_bytes(source),
        "replayed_scores_sha256": sha256_bytes(replayed),
        "byte_identical_replay": replayed == source,
        "passed": all(record_checks) and replayed == source,
        "mpra_functional_outcomes_accessed": False,
    }
    args.replay_receipt_json.write_text(
        json.dumps(replay_receipt, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"cache_records={len(cache_entries)} "
        f"bytes={cache_manifest['total_bytes']} "
        f"replay={'pass' if replay_receipt['passed'] else 'fail'}"
    )
    if not replay_receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
