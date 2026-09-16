#!/usr/bin/env python3
"""Take the first 1,000 uniquely lifted candidates by frozen selection rank."""

from __future__ import annotations

import argparse
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
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def set_sha256(rows: list[dict[str, str]], id_field: str) -> str:
    payload = "\n".join(sorted(row[id_field] for row in rows)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_hg19_tsv", type=Path)
    parser.add_argument("lifted_hg38_tsv", type=Path)
    parser.add_argument("final_hg19_tsv", type=Path)
    parser.add_argument("final_hg38_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    parser.add_argument("--target", type=int, default=1000)
    parser.add_argument("--reference-validation-tsv", type=Path)
    args = parser.parse_args()

    candidates = read_tsv(args.candidate_hg19_tsv)
    lifted_all = read_tsv(args.lifted_hg38_tsv)
    lifted = lifted_all
    valid_hg38_ids = None
    if args.reference_validation_tsv is not None:
        validation = read_tsv(args.reference_validation_tsv)
        valid_hg38_ids = {
            row["variant_id"]
            for row in validation
            if row["reference_status"] == "match"
        }
        lifted = [row for row in lifted if row["variant_id"] in valid_hg38_ids]
    lifted_by_source = {row["source_variant_id_hg19"]: row for row in lifted}
    ranked = sorted(candidates, key=lambda row: int(row["selection_rank"]))
    selected_hg19 = [
        row for row in ranked if row["variant_id"] in lifted_by_source
    ][: args.target]
    if len(selected_hg19) != args.target:
        raise ValueError(
            f"only {len(selected_hg19)} uniquely lifted candidates; "
            f"expected {args.target}"
        )
    selected_hg38 = [lifted_by_source[row["variant_id"]] for row in selected_hg19]
    selected_ids = {row["variant_id"] for row in selected_hg19}
    failed_before_target = [
        row["variant_id"]
        for row in ranked
        if int(row["selection_rank"]) <= args.target
        and row["variant_id"] not in lifted_by_source
    ]
    reserves_used = sum(
        int(row["selection_rank"]) > args.target for row in selected_hg19
    )

    write_tsv(args.final_hg19_tsv, selected_hg19)
    write_tsv(args.final_hg38_tsv, selected_hg38)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_variants": len(candidates),
        "uniquely_lifted_candidates": len(lifted_all),
        "reference_valid_lifted_candidates": len(lifted),
        "reference_validation_applied": args.reference_validation_tsv is not None,
        "final_variants": len(selected_hg19),
        "failed_primary_candidates": failed_before_target,
        "reserves_used": reserves_used,
        "first_excluded_valid_reserve_rank": min(
            int(row["selection_rank"])
            for row in lifted
            if row["source_variant_id_hg19"] not in selected_ids
        ),
        "final_hg19_set_sha256": set_sha256(selected_hg19, "variant_id"),
        "final_hg38_set_sha256": set_sha256(selected_hg38, "variant_id"),
        "mpra_functional_outcomes_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"final={len(selected_hg19)} reserves_used={reserves_used} "
        f"failed_primary={len(failed_before_target)}"
    )


if __name__ == "__main__":
    main()
