#!/usr/bin/env python3
"""Finalize the pilot while preserving every frozen stratum target."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def set_hash(rows: list[dict[str, str]], field: str) -> str:
    payload = "\n".join(sorted(row[field] for row in rows)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_hg19_tsv", type=Path)
    parser.add_argument("lifted_hg38_tsv", type=Path)
    parser.add_argument("reference_validation_tsv", type=Path)
    parser.add_argument("final_hg19_tsv", type=Path)
    parser.add_argument("final_hg38_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    args = parser.parse_args()

    candidates = read_tsv(args.candidate_hg19_tsv)
    lifted = read_tsv(args.lifted_hg38_tsv)
    validation = read_tsv(args.reference_validation_tsv)
    valid_hg38 = {row["variant_id"] for row in validation if row["reference_status"] == "match"}
    lifted_by_source = {
        row["source_variant_id_hg19"]: row for row in lifted if row["variant_id"] in valid_hg38
    }

    by_stratum: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for row in candidates:
        by_stratum[row["stratum_id"]].append(row)

    final_hg19: list[dict[str, str]] = []
    shortfalls = []
    for stratum_id, rows in sorted(by_stratum.items()):
        ranked = sorted(rows, key=lambda row: int(row["stratum_rank"]))
        target = int(rows[0]["stratum_target"])
        selected = [row for row in ranked if row["variant_id"] in lifted_by_source][:target]
        if len(selected) != target:
            shortfalls.append({"stratum_id": stratum_id, "target": target, "available": len(selected)})
        final_hg19.extend(selected)
    if shortfalls:
        raise ValueError(f"frozen stratum shortfalls: {shortfalls}")
    if len(final_hg19) != 10_000:
        raise ValueError(f"final variants={len(final_hg19)}, expected 10000")

    final_hg19.sort(key=lambda row: (int(row["selection_rank"]), row["variant_id"]))
    final_hg38 = [lifted_by_source[row["variant_id"]] for row in final_hg19]
    selected_ids = {row["variant_id"] for row in final_hg19}
    failed_primary = [
        row["variant_id"] for row in candidates
        if row["selection_tier"] == "primary" and row["variant_id"] not in lifted_by_source
    ]
    reserves_used = [
        row["variant_id"] for row in final_hg19 if row["selection_tier"] == "reserve"
    ]

    write_tsv(args.final_hg19_tsv, final_hg19)
    write_tsv(args.final_hg38_tsv, final_hg38)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_variants": len(candidates),
        "uniquely_lifted_candidates": len(lifted),
        "reference_valid_lifted_candidates": len(lifted_by_source),
        "final_variants": len(final_hg19),
        "failed_primary_candidates": failed_primary,
        "reserves_used": len(reserves_used),
        "reserve_variant_ids": reserves_used,
        "stratum_shortfalls": shortfalls,
        "frozen_stratum_targets_preserved": not shortfalls,
        "final_hg19_set_sha256": set_hash(final_hg19, "variant_id"),
        "final_hg38_set_sha256": set_hash(final_hg38, "variant_id"),
        "all_final_hg38_reference_valid": all(row["variant_id"] in valid_hg38 for row in final_hg38),
        "mpra_functional_outcomes_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"final=10000 failed_primary={len(failed_primary)} reserves_used={len(reserves_used)}")


if __name__ == "__main__":
    main()
