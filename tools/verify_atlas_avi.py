#!/usr/bin/env python3
"""Independently inspect outcome-blind AVI baseline artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("avi_dir", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--min-complete-fraction", type=float, default=0.99)
    args = parser.parse_args()
    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    with (args.avi_dir / "avi_status.tsv").open(newline="") as handle:
        status = list(csv.DictReader(handle, delimiter="\t"))
    scores = [json.loads(line) for line in (args.avi_dir / "avi_scores.jsonl").read_text().splitlines()]
    summary = json.loads((args.avi_dir / "avi_summary.json").read_text())
    manifest_ids = {row["variant_id"] for row in manifest}
    status_ids = {row["variant_id_hg38"] for row in status}
    complete_ids = {row["variant_id_hg38"] for row in status if row["status"] == "complete"}
    score_ids = {row["variant_id_hg38"] for row in scores}
    complete_fraction = len(complete_ids) / len(manifest_ids)
    checks = {
        "manifest_has_10000_unique_variants": len(manifest) == len(manifest_ids) == 10_000,
        "status_ids_equal_manifest": status_ids == manifest_ids and len(status) == 10_000,
        "score_ids_equal_complete_ids": score_ids == complete_ids and len(scores) == len(score_ids),
        "scores_are_finite": all(
            math.isfinite(float(row[field])) for row in scores for field in ("raw_score", "quantile_score")
        ),
        "completion_meets_minimum": complete_fraction >= args.min_complete_fraction,
        "summary_completion_matches": math.isclose(summary["complete_fraction"], complete_fraction),
        "summary_checksum_matches": summary["scores_sha256"] == sha256(args.avi_dir / "avi_scores.jsonl"),
        "unsigned_avi_scorer": summary["scorer_is_signed"] is False,
        "no_mpra_outcomes_accessed": summary["mpra_functional_outcomes_accessed"] is False,
        "no_rule_mutation": summary["scientific_rules_mutated"] is False,
    }
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "passed": all(checks.values()),
        "complete_fraction": complete_fraction,
        "complete_variants": len(complete_ids),
        "scores_sha256": sha256(args.avi_dir / "avi_scores.jsonl"),
        "mpra_functional_outcomes_accessed": False,
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"avi_inspection={'pass' if receipt['passed'] else 'fail'} complete={len(complete_ids)}")
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
