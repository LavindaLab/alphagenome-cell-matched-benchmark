#!/usr/bin/env python3
"""Consolidate the outcome-blind engineering gate and scale projections."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FULL_CEILING = 281_327
LOCKED_PILOT = 10_000


def load_json(path: Path):
    return json.loads(path.read_text())


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("engineering_dir", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("projection_json", type=Path)
    args = parser.parse_args()
    root = args.engineering_dir
    run1 = root / "run1"
    run2 = root / "run2"

    freeze = load_json(root / "PREQUERY_FREEZE.json")
    finalization = load_json(root / "finalization_manifest.json")
    final_variants = read_tsv(root / "variants_hg38.tsv")
    validation = read_tsv(root / "candidate_reference_validation_hg38.tsv")
    valid_ids = {
        row["variant_id"]
        for row in validation
        if row["reference_status"] == "match"
    }
    final_ids = {row["variant_id"] for row in final_variants}
    summary1 = load_json(run1 / "atlas_retrieval_summary.json")
    summary2 = load_json(run2 / "atlas_retrieval_summary.json")
    inspection1 = load_json(run1 / "inspection_receipt.json")
    inspection2 = load_json(run2 / "inspection_receipt.json")
    cache = load_json(root / "cache_manifest.json")
    cache_replay = load_json(root / "cache_replay_receipt.json")
    features = load_json(root / "prediction_features_summary.json")
    selection = load_json(root / "candidate_selection_manifest.json")

    protocol_path = root.parent / "PROTOCOL_v0.1.md"
    ledger_path = root.parent / "EVIDENCE_LEDGER.md"
    scores1 = run1 / "atlas_scores.jsonl"
    scores2 = run2 / "atlas_scores.jsonl"
    checks = {
        "prequery_selection_was_outcome_blind": freeze[
            "mpra_functional_outcomes_accessed"
        ]
        is False,
        "protocol_hash_unchanged": sha256(protocol_path)
        == freeze["protocol_sha256"],
        "evidence_ledger_hash_unchanged": sha256(ledger_path)
        == freeze["evidence_ledger_sha256"],
        "final_manifest_has_1000_unique_variants": (
            len(final_variants) == len(final_ids) == 1000
        ),
        "all_final_variants_are_reference_valid": final_ids <= valid_ids,
        "coordinate_failures_replaced_only_by_frozen_reserve": (
            finalization["reserves_used"]
            == len(finalization["failed_primary_candidates"])
            == 4
        ),
        "run1_independent_inspection_passed": inspection1["passed"],
        "run2_independent_inspection_passed": inspection2["passed"],
        "both_runs_100_percent_complete": (
            summary1["complete_variant_fraction"]
            == summary2["complete_variant_fraction"]
            == 1.0
        ),
        "both_runs_used_zero_retries": (
            summary1["variant_request_attempts"]
            == summary2["variant_request_attempts"]
            == 1000
        ),
        "normalized_full_rerun_is_byte_identical": (
            scores1.read_bytes() == scores2.read_bytes()
            and summary1["normalized_scores_sha256"]
            == summary2["normalized_scores_sha256"]
        ),
        "cache_replay_is_byte_identical": cache_replay["passed"],
        "cache_has_one_record_per_variant": cache["records"] == 1000,
        "prediction_plumbing_has_5000_variant_cell_rows": (
            features["variant_cell_rows"] == 5000
            and features["unique_variants"] == 1000
        ),
        "no_coverage_threshold_or_ranking_was_applied": (
            features["coverage_threshold_applied"] is False
            and features["ranking_selection_applied"] is False
        ),
        "no_mpra_outcomes_were_accessed": all(
            item["mpra_functional_outcomes_accessed"] is False
            for item in (summary1, summary2, cache, cache_replay, features)
        ),
        "no_scientific_rules_were_mutated": (
            summary1["scientific_rules_mutated"] is False
            and summary2["scientific_rules_mutated"] is False
            and features["scientific_rules_mutated"] is False
        ),
    }

    mean_seconds = (summary1["elapsed_seconds"] + summary2["elapsed_seconds"]) / 2
    scalable_bytes_1000 = (
        (run1 / "atlas_scores.jsonl").stat().st_size
        + (run1 / "atlas_variant_status.tsv").stat().st_size
    )
    fixed_run_bytes = sum(
        path.stat().st_size for path in run1.iterdir() if path.is_file()
    ) - scalable_bytes_1000
    cache_bytes_1000 = cache["total_bytes"]
    feature_bytes_1000 = (root / "prediction_features.tsv").stat().st_size

    def projection(variants: int):
        linear_seconds = mean_seconds * variants / 1000
        normalized_bytes = scalable_bytes_1000 * variants / 1000 + fixed_run_bytes
        cache_bytes = cache_bytes_1000 * variants / 1000
        feature_bytes = feature_bytes_1000 * variants / 1000
        return {
            "unique_variant_requests": variants,
            "metadata_requests": 1,
            "total_api_requests_without_retries": variants + 1,
            "linear_wall_seconds_at_5_workers": round(linear_seconds, 2),
            "two_to_three_x_wall_sensitivity_seconds": [
                round(linear_seconds * 2, 2),
                round(linear_seconds * 3, 2),
            ],
            "normalized_scores_status_and_fixed_metadata_bytes": round(
                normalized_bytes
            ),
            "per_variant_cache_bytes": round(cache_bytes),
            "prediction_feature_table_bytes": round(feature_bytes),
            "combined_normalized_cache_and_feature_bytes": round(
                normalized_bytes + cache_bytes + feature_bytes
            ),
        }

    projection_output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "basis": {
            "workers": 5,
            "run1_seconds": summary1["elapsed_seconds"],
            "run2_seconds": summary2["elapsed_seconds"],
            "mean_seconds_per_1000": round(mean_seconds, 3),
            "observed_variant_requests_per_second": round(1000 / mean_seconds, 3),
            "scalable_normalized_bytes_per_1000": scalable_bytes_1000,
            "fixed_run_metadata_bytes": fixed_run_bytes,
            "cache_bytes_per_1000": cache_bytes_1000,
            "prediction_feature_bytes_per_1000": feature_bytes_1000,
            "retry_rate": 0.0,
        },
        "locked_10000_pilot": projection(LOCKED_PILOT),
        "confirmation_split_all_snv_ceiling": projection(
            selection["eligible_unique_canonical_snvs_all_splits"]
            - selection["eligible_unique_canonical_snvs_gate_blocks"]
        ),
        "full_281327_snv_ceiling": projection(FULL_CEILING),
        "limitations": [
            "linear timing assumes the observed five-worker rate persists",
            "two-to-three-x sensitivity is operational headroom, not a measured SLA",
            "storage projects normalized outputs and cache records, not raw protobufs",
            "full ceiling precedes noncoding and other frozen eligibility exclusions",
            "one rerun would approximately double API time and normalized run storage",
        ],
        "mpra_functional_outcomes_accessed": False,
    }
    args.projection_json.write_text(
        json.dumps(projection_output, indent=2, sort_keys=True) + "\n"
    )

    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "outcome-blind 1,000-variant engineering/scaling gate",
        "checks": checks,
        "passed": all(checks.values()),
        "run1_scores_sha256": sha256(scores1),
        "run2_scores_sha256": sha256(scores2),
        "projected_requirements_path": args.projection_json.name,
        "mpra_functional_outcomes_accessed": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        f"engineering_inspection={'pass' if receipt['passed'] else 'fail'} "
        f"requests_per_second={projection_output['basis']['observed_variant_requests_per_second']}"
    )
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
