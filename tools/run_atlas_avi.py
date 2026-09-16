#!/usr/bin/env python3
"""Retrieve the prespecified unsigned AVI baseline without MPRA outcome access."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from run_atlas_metadata_check import file_sha256, import_sdk, json_safe, load_variants, metadata_records


SCORER = "AVI_SCORE"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variants_hg38_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--sdk-path", type=Path)
    parser.add_argument("--expected-variants", type=int, default=10_000)
    parser.add_argument("--max-workers", type=int, default=5)
    parser.add_argument("--retries", type=int, default=2)
    args = parser.parse_args()
    rows = load_variants(args.variants_hg38_tsv, args.expected_variants)
    atlas, genome, _, sdk_version = import_sdk(args.sdk_path)
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise SystemExit("ALPHA_GENOME_API_KEY is not set")
    client = atlas.create(api_key, timeout=30)
    metadata = client.scorer_metadata()
    if SCORER not in metadata:
        raise RuntimeError(f"required scorer absent: {SCORER}")
    scorer_info = metadata[SCORER]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    preflight = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "prespecified outcome-blind AVI baseline retrieval",
        "variant_manifest_sha256": file_sha256(args.variants_hg38_tsv),
        "variants": len(rows),
        "requested_scorer": SCORER,
        "scorer_is_signed": scorer_info.is_signed,
        "track_metadata": metadata_records(scorer_info.track_metadata),
        "sdk_version": sdk_version,
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    (args.output_dir / "avi_preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    started = time.monotonic()

    def retrieve(row: dict[str, str]):
        variant = genome.Variant(
            chromosome=row["chromosome"],
            position=int(row["position_1based"]),
            reference_bases=row["reference"],
            alternate_bases=row["alternate"],
        )
        error = ""
        for attempt in range(1, args.retries + 2):
            try:
                response = client.query_variant(variant, requested_scorers=[SCORER])
                data = response[SCORER]
                quantiles = data.layers.get("quantiles")
                raw = json_safe(data.X[0, 0])
                quantile = json_safe(quantiles[0, 0]) if quantiles is not None else None
                valid = all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in (raw, quantile)
                )
                result = {
                    "variant_id_hg38": row["variant_id"],
                    "source_variant_id_hg19": row.get("source_variant_id_hg19", ""),
                    "scorer_name": SCORER,
                    "raw_score": raw,
                    "quantile_score": quantile,
                }
                return row, result if valid else None, attempt, "." if valid else "non-finite score"
            except Exception as exc:
                error = " ".join(f"{type(exc).__name__}: {exc}".split())
                if attempt <= args.retries:
                    time.sleep(min(2 ** (attempt - 1), 4))
        return row, None, args.retries + 1, error

    retrievals = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [executor.submit(retrieve, row) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            retrievals.append(future.result())
    scores = sorted(
        (result for _, result, _, _ in retrievals if result is not None),
        key=lambda row: row["variant_id_hg38"],
    )
    status = sorted(
        (
            {
                "variant_id_hg38": row["variant_id"],
                "source_variant_id_hg19": row.get("source_variant_id_hg19", ""),
                "status": "complete" if result is not None else "incomplete",
                "attempts": attempts,
                "error": error,
            }
            for row, result, attempts, error in retrievals
        ),
        key=lambda row: row["variant_id_hg38"],
    )
    lines = [json.dumps(row, sort_keys=True) for row in scores]
    score_bytes = ("\n".join(lines) + "\n").encode()
    (args.output_dir / "avi_scores.jsonl").write_bytes(score_bytes)
    with (args.output_dir / "avi_status.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(status[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(status)
    counts = Counter(row["status"] for row in status)
    summary = {
        **preflight,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "status_counts": dict(sorted(counts.items())),
        "complete_fraction": counts["complete"] / len(status),
        "variant_request_attempts": sum(int(row["attempts"]) for row in status),
        "scores": len(scores),
        "scores_sha256": hashlib.sha256(score_bytes).hexdigest(),
        "api_key_persisted": False,
    }
    (args.output_dir / "avi_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"variants={len(rows)} complete={counts['complete']} incomplete={counts['incomplete']} requests={summary['variant_request_attempts']}")


if __name__ == "__main__":
    main()
