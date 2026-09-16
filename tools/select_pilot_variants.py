#!/usr/bin/env python3
"""Freeze the outcome-blind, proportionally stratified 10,000-SNV pilot."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from inventory_pilot_strata import (
    BLOCK_SALT,
    sha256,
    table1_metadata,
    table4_masks,
)


TARGET = 10_000
SELECTION_SALT = "AGATLAS-PILOT-10000-2026-V01"


def digest(label: str) -> str:
    return hashlib.sha256(f"{SELECTION_SALT}\0{label}".encode()).hexdigest()


def block_fields(variant: str) -> tuple[str, str]:
    chromosome, position, _, _ = variant.split(":")
    block_id = f"{chromosome}:{(int(position) - 1) // 1_000_000}"
    block_hash = hashlib.sha256(f"{BLOCK_SALT}\0{block_id}".encode()).hexdigest()
    return block_id, block_hash


def allocate(strata: dict[str, list[str]], target: int) -> dict[str, int]:
    total = sum(map(len, strata.values()))
    quotas = {key: target * len(values) / total for key, values in strata.items()}
    result = {key: math.floor(value) for key, value in quotas.items()}
    remaining = target - sum(result.values())
    order = sorted(
        strata,
        key=lambda key: (-(quotas[key] - result[key]), digest(f"stratum:{key}"), key),
    )
    for key in order[:remaining]:
        result[key] += 1
    assert sum(result.values()) == target
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table4", type=Path)
    parser.add_argument("table1", type=Path)
    parser.add_argument("candidates_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    parser.add_argument(
        "--full-reserve-stratum",
        action="append",
        default=[],
        help="pre-outcome technical amendment: freeze every ranked reserve in this stratum",
    )
    args = parser.parse_args()

    measurement_masks = table4_masks(args.table4)
    categories, _ = table1_metadata(args.table1, set(measurement_masks))

    metadata: dict[str, dict[str, str]] = {}
    strata: dict[str, list[str]] = collections.defaultdict(list)
    for variant, mask in measurement_masks.items():
        category_set = "|".join(sorted(categories.get(variant, {"NOT_IN_TABLE1"})))
        chromosome = variant.split(":", 1)[0]
        stratum_id = f"{chromosome}::{mask}::{category_set}"
        metadata[variant] = {
            "measurement_mask": mask,
            "siraj_category_set": category_set,
            "stratum_id": stratum_id,
        }
        strata[stratum_id].append(variant)

    targets = allocate(strata, TARGET)
    selected: list[dict[str, str | int]] = []
    allocations = []
    for stratum_id, variants in sorted(strata.items()):
        ranked = sorted(variants, key=lambda variant: (digest(variant), variant))
        target = targets[stratum_id]
        if stratum_id in args.full_reserve_stratum:
            reserve_count = len(ranked) - target
        else:
            reserve_count = min(
                len(ranked) - target, max(2, math.ceil(target * 0.05))
            )
        candidate_count = target + reserve_count
        for rank, variant in enumerate(ranked[:candidate_count], start=1):
            chromosome, position, reference, alternate = variant.split(":")
            block_id, block_hash = block_fields(variant)
            selected.append(
                {
                    "variant_id": variant,
                    "chromosome": chromosome,
                    "position_1based": int(position),
                    "reference": reference,
                    "alternate": alternate,
                    "block_id_1mb": block_id,
                    "block_sha256": block_hash,
                    "frozen_split": "gate",
                    **metadata[variant],
                    "stratum_population": len(ranked),
                    "stratum_target": target,
                    "stratum_rank": rank,
                    "selection_tier": "primary" if rank <= target else "reserve",
                    "selection_sha256": digest(variant),
                }
            )
        allocations.append(
            {
                "stratum_id": stratum_id,
                "population": len(ranked),
                "target": target,
                "candidate_reserves": reserve_count,
            }
        )

    selected.sort(
        key=lambda row: (
            row["selection_tier"] != "primary",
            row["selection_sha256"],
            row["variant_id"],
        )
    )
    for rank, row in enumerate(selected, start=1):
        row["selection_rank"] = rank

    args.candidates_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.candidates_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(selected[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected)

    sample_payload = "\n".join(sorted(row["variant_id"] for row in selected)) + "\n"
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind frozen 10,000-SNV scientific pilot selection",
        "protocol": "PROTOCOL_v0.1.md",
        "table4_sha256": sha256(args.table4),
        "table1_sha256": sha256(args.table1),
        "block_salt": BLOCK_SALT,
        "selection_salt": SELECTION_SALT,
        "gate_rule": "integer SHA-256 block digest is in lowest 10% of hash space",
        "stratum_definition": "chromosome x cell-measurement-availability mask x exact set of Siraj Category values",
        "allocation": "proportional largest remainder; remainder ties by salted stratum SHA-256",
        "within_stratum_ranking": "salted variant SHA-256 then normalized variant identifier",
        "reserve_rule": "per stratum, max(2, ceiling(5% of target)), capped by remaining population; every remaining frozen-ranked candidate is included for strata named in full_reserve_strata",
        "full_reserve_strata": sorted(args.full_reserve_stratum),
        "gate_canonical_snv_population": len(measurement_masks),
        "strata": len(strata),
        "target_variants": TARGET,
        "primary_candidates": sum(row["selection_tier"] == "primary" for row in selected),
        "reserve_candidates": sum(row["selection_tier"] == "reserve" for row in selected),
        "candidate_variants": len(selected),
        "candidate_sample_sha256": hashlib.sha256(sample_payload.encode()).hexdigest(),
        "stratum_allocations": allocations,
        "table4_fields_used": ["Variant", "cell-specific activity-estimate missingness only"],
        "table1_fields_used": ["Variant", "Category"],
        "numeric_mpra_values_retained": False,
        "mpra_activity_calls_accessed": False,
        "mpra_emvar_calls_accessed": False,
        "mpra_skew_values_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"population={len(measurement_masks)} strata={len(strata)} "
        f"primary={manifest['primary_candidates']} reserves={manifest['reserve_candidates']}"
    )


if __name__ == "__main__":
    main()
