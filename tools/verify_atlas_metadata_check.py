#!/usr/bin/env python3
"""Independently verify Atlas feasibility outputs without MPRA outcome access."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


CONTEXTS = {
    "EFO:0001086",
    "EFO:0001187",
    "EFO:0002067",
    "EFO:0002824",
    "EFO:0003072",
}
MODALITIES = {"DNASE", "H3K27ac"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("--expected-variants", type=int, default=100)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--min-complete-fraction", type=float, default=1.0)
    args = parser.parse_args()
    root = args.artifact_dir

    manifest_path = args.manifest or (root / "variants_hg38.tsv")
    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    manifest_ids = {row["variant_id"] for row in manifest}

    score_path = root / "atlas_scores.jsonl"
    scores = [json.loads(line) for line in score_path.read_text().splitlines()]
    scores_by_variant: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in scores:
        scores_by_variant[str(row["variant_id_hg38"])].append(row)

    expected_pairs = set(itertools.product(CONTEXTS, MODALITIES))
    exact_pairs = {}
    finite_scores = {}
    for variant, rows in scores_by_variant.items():
        pairs = [(str(row["ontology_curie"]), str(row["modality"])) for row in rows]
        exact_pairs[variant] = len(pairs) == 10 and set(pairs) == expected_pairs
        finite_scores[variant] = all(
            isinstance(row[field], (int, float))
            and not isinstance(row[field], bool)
            and math.isfinite(float(row[field]))
            for row in rows
            for field in ("raw_score", "quantile_score")
        )

    with (root / "atlas_variant_status.tsv").open(newline="") as handle:
        status_rows = list(csv.DictReader(handle, delimiter="\t"))
    status_counts = Counter(row["status"] for row in status_rows)
    status_ids = {row["variant_id_hg38"] for row in status_rows}
    complete_ids = {
        row["variant_id_hg38"] for row in status_rows if row["status"] == "complete"
    }
    incomplete_ids = status_ids - complete_ids
    complete_fraction = len(complete_ids) / args.expected_variants
    summary = json.loads((root / "atlas_retrieval_summary.json").read_text())
    metadata = json.loads((root / "atlas_scorer_metadata.json").read_text())

    checks = {
        "manifest_has_expected_unique_variants": (
            len(manifest) == len(manifest_ids) == args.expected_variants
        ),
        "status_variants_equal_manifest": status_ids == manifest_ids,
        "score_variants_equal_complete_status_variants": set(scores_by_variant)
        == complete_ids,
        "expected_complete_score_row_count": len(scores) == len(complete_ids) * 10,
        "each_complete_variant_has_exactly_one_row_per_required_context_modality": all(
            exact_pairs.values()
        ),
        "all_raw_and_quantile_scores_are_finite": all(finite_scores.values()),
        "all_status_rows_are_complete_or_incomplete": len(status_rows)
        == args.expected_variants
        and set(status_counts) <= {"complete", "incomplete"},
        "incomplete_variants_have_no_score_rows": not (
            incomplete_ids & set(scores_by_variant)
        ),
        "completion_meets_prespecified_minimum": complete_fraction
        >= args.min_complete_fraction,
        "summary_completion_fraction_matches": math.isclose(
            summary.get("complete_variant_fraction", -1), complete_fraction
        ),
        "summary_score_checksum_matches_file": summary.get(
            "normalized_scores_sha256"
        )
        == sha256(score_path),
        "both_scorers_have_all_five_contexts": set(metadata) == MODALITIES
        and all(item.get("all_contexts_present") for item in metadata.values()),
        "no_mpra_outcomes_accessed": summary.get(
            "mpra_functional_outcomes_accessed"
        )
        is False,
        "no_scientific_rules_mutated": summary.get("scientific_rules_mutated")
        is False,
    }
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "inspection_scope": "Atlas schema/context/retrieval artifacts only; no MPRA workbook or outcomes",
        "checks": checks,
        "passed": all(checks.values()),
        "manifest_variants": len(manifest),
        "score_rows": len(scores),
        "status_counts": dict(status_counts),
        "atlas_scores_sha256": sha256(score_path),
        "atlas_summary_sha256": sha256(root / "atlas_retrieval_summary.json"),
        "mpra_functional_outcomes_accessed": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        f"inspection={'pass' if receipt['passed'] else 'fail'} "
        f"variants={len(manifest)} score_rows={len(scores)}"
    )
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
