#!/usr/bin/env python3
"""Verify that outcome extraction is complete and confined to frozen gate IDs."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


CELLS = {"A549", "HepG2", "K562", "SK-N-SH", "HCT116"}


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("final_manifest", type=Path)
    parser.add_argument("outcomes", type=Path)
    parser.add_argument("variant_metadata", type=Path)
    parser.add_argument("noncoding", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    final = rows(args.final_manifest)
    outcomes = rows(args.outcomes)
    metadata = rows(args.variant_metadata)
    noncoding = rows(args.noncoding)
    final_hg38 = {row["variant_id"] for row in final}
    final_hg19 = {row["source_variant_id_hg19"] for row in final}
    outcome_pairs = {(row["variant_id_hg38"], row["cell_context"]) for row in outcomes}
    bool_values = {"", "NA", "true", "false"}
    checks = {
        "final_manifest_is_10000_gate_variants": len(final) == len(final_hg38) == len(final_hg19) == 10_000 and all(row["frozen_split"] == "gate" for row in final),
        "outcomes_are_exact_10000_by_5_cartesian_product": len(outcomes) == len(outcome_pairs) == 50_000 and outcome_pairs == {(variant, cell) for variant in final_hg38 for cell in CELLS},
        "outcome_source_ids_equal_gate_manifest": {row["source_variant_id_hg19"] for row in outcomes} == final_hg19,
        "outcome_boolean_schema_valid": all(row[field] in bool_values for row in outcomes for field in ("measured", "active", "enhancer_active", "emvar")),
        "enhancer_active_implies_active_and_positive_activity": all(
            row["enhancer_active"] != "true" or (row["active"] == "true" and float(row["activity_log2fc"]) > 0)
            for row in outcomes
        ),
        "emvar_implies_measured": all(row["emvar"] != "true" or row["measured"] == "true" for row in outcomes),
        "variant_metadata_exactly_matches_gate": len(metadata) == 10_000 and {row["variant_id_hg38"] for row in metadata} == final_hg38,
        "noncoding_exactly_matches_gate": len(noncoding) == 10_000 and {row["variant_id_hg38"] for row in noncoding} == final_hg38,
    }
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "passed": all(checks.values()),
        "outcome_rows_by_cell": dict(sorted(Counter(row["cell_context"] for row in outcomes).items())),
        "confirmation_split_accessed": False,
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"gate_outcome_inspection={'pass' if receipt['passed'] else 'fail'} checks={len(checks)}")
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
