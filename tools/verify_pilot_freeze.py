#!/usr/bin/env python3
"""Independently verify the outcome-blind 10,000-SNV pilot freeze."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ALLOWED_COLUMNS = {
    "variant_id", "chromosome", "position_1based", "reference", "alternate",
    "block_id_1mb", "block_sha256", "frozen_split", "measurement_mask",
    "siraj_category_set", "stratum_id", "stratum_population", "stratum_target",
    "stratum_rank", "selection_tier", "selection_sha256", "selection_rank",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pilot_dir", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("--candidate-file", default="candidate_variants_hg19.tsv")
    parser.add_argument("--selection-file", default="selection_manifest.json")
    parser.add_argument("--freeze-file", default="PRE_OUTCOME_FREEZE.json")
    args = parser.parse_args()
    root = args.pilot_dir
    candidate_path = root / args.candidate_file
    selection_path = root / args.selection_file
    freeze_path = root / args.freeze_file
    with candidate_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = json.loads(selection_path.read_text())
    freeze = json.loads(freeze_path.read_text())

    by_stratum: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for row in rows:
        by_stratum[row["stratum_id"]].append(row)
    allocation = {item["stratum_id"]: item for item in manifest["stratum_allocations"]}
    checks = {
        "candidate_manifest_hash_frozen": sha256(candidate_path) == freeze["candidate_manifest_sha256"],
        "selection_manifest_hash_frozen": sha256(selection_path) == freeze["selection_manifest_sha256"],
        "candidate_ids_unique": len(rows) == len({row["variant_id"] for row in rows}) == manifest["candidate_variants"],
        "only_allowed_columns": set(rows[0]) == ALLOWED_COLUMNS,
        "exactly_10000_primary_candidates": sum(row["selection_tier"] == "primary" for row in rows) == 10_000,
        "all_candidates_are_gate_split": all(row["frozen_split"] == "gate" for row in rows),
        "stratum_count_matches": len(by_stratum) == manifest["strata"] == 362,
        "stratum_primary_targets_match": all(
            sum(row["selection_tier"] == "primary" for row in stratum_rows) == allocation[stratum_id]["target"]
            for stratum_id, stratum_rows in by_stratum.items()
        ),
        "no_outcome_access_recorded": all(
            freeze[key] is False for key in (
                "mpra_numeric_values_retained", "mpra_activity_calls_accessed",
                "mpra_emvar_calls_accessed", "mpra_skew_values_accessed",
                "atlas_pilot_predictions_accessed", "confirmation_data_accessed",
            )
        ),
        "no_scientific_rule_mutation": freeze["scientific_rules_mutated"] is False,
    }
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "pre-outcome 10,000-SNV pilot selection freeze",
        "checks": checks,
        "passed": all(checks.values()),
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"pilot_freeze_inspection={'pass' if receipt['passed'] else 'fail'} checks={len(checks)}")
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
