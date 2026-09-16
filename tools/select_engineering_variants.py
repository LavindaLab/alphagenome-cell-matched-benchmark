#!/usr/bin/env python3
"""Freeze the outcome-blind 1,000-SNV engineering sample from gate blocks."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from select_metadata_variants import (
    CHROMOSOMES,
    PRE_OUTCOME_COMMIT,
    extract_variant_strings,
    file_sha256,
    sample_sha256,
    selection_digest,
)


BLOCK_SALT = "AGATLAS-SIRAJ-2026-V01"
SELECTION_SALT = "AGATLAS-ENGINEERING-1000-2026-V01"
TARGET = 1000
CANDIDATES = 1020
PER_CHROMOSOME = 40


def block_fields(variant: str) -> tuple[str, str, bool]:
    chromosome, position, _, _ = variant.split(":")
    block_id = f"{chromosome}:{(int(position) - 1) // 1_000_000}"
    digest = hashlib.sha256(f"{BLOCK_SALT}\0{block_id}".encode()).hexdigest()
    in_gate = int(digest, 16) * 10 < 2**256
    return block_id, digest, in_gate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_xlsx", type=Path)
    parser.add_argument("variants_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    args = parser.parse_args()

    all_variants = extract_variant_strings(args.source_xlsx)
    gate_variants = {
        variant for variant in all_variants if block_fields(variant)[2]
    }
    by_chromosome: dict[str, list[str]] = collections.defaultdict(list)
    for variant in gate_variants:
        by_chromosome[variant.split(":", 1)[0]].append(variant)

    selected: list[tuple[str, str]] = []
    for chromosome in CHROMOSOMES:
        ranked = sorted(
            by_chromosome[chromosome],
            key=lambda item: (selection_digest(SELECTION_SALT, item), item),
        )
        selected.extend(
            (variant, "chromosome_quota")
            for variant in ranked[:PER_CHROMOSOME]
        )
    selected_ids = {variant for variant, _ in selected}
    remaining = sorted(
        gate_variants - selected_ids,
        key=lambda item: (selection_digest(SELECTION_SALT, item), item),
    )
    selected.extend(
        (variant, "global_fill")
        for variant in remaining[: CANDIDATES - len(selected)]
    )
    if len(selected) != CANDIDATES:
        raise ValueError(
            f"selected {len(selected)} candidates; expected {CANDIDATES}"
        )

    rows = []
    rank_by_variant = {
        variant: rank for rank, (variant, _) in enumerate(selected, start=1)
    }
    for variant, tier in sorted(
        selected,
        key=lambda pair: (
            CHROMOSOMES.index(pair[0].split(":", 1)[0]),
            int(pair[0].split(":")[1]),
            pair[0],
        ),
    ):
        chromosome, position, reference, alternate = variant.split(":")
        block_id, block_sha256, in_gate = block_fields(variant)
        assert in_gate
        rows.append(
            {
                "variant_id": variant,
                "chromosome": chromosome,
                "position_1based": int(position),
                "reference": reference,
                "alternate": alternate,
                "block_id_1mb": block_id,
                "block_sha256": block_sha256,
                "frozen_split": "gate",
                "selection_tier": tier,
                "selection_rank": rank_by_variant[variant],
                "selection_sha256": selection_digest(SELECTION_SALT, variant),
            }
        )

    args.variants_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.variants_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    chromosome_counts = collections.Counter(row["chromosome"] for row in rows)
    block_counts = collections.Counter(row["block_id_1mb"] for row in rows)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind 1,000-variant Atlas engineering/scaling check",
        "pre_outcome_commit": PRE_OUTCOME_COMMIT,
        "source_name": args.source_xlsx.name,
        "source_sha256": file_sha256(args.source_xlsx),
        "xlsx_members_read_for_selection": ["xl/sharedStrings.xml"],
        "worksheet_cells_read_for_selection": False,
        "mpra_functional_outcomes_accessed_for_selection": False,
        "block_definition": "chromosome:floor((position-1)/1,000,000)",
        "block_hash_algorithm": "SHA-256(block_salt + NUL + block_id)",
        "block_salt": BLOCK_SALT,
        "gate_rule": "integer SHA-256 digest is in the lowest 10% of hash space",
        "selection_salt": SELECTION_SALT,
        "selection_algorithm": (
            "within gate blocks only, up to 40 lowest salted SHA-256 canonical "
            "SNVs per chr1-22/chrX, then globally lowest unused SNVs; freeze "
            "1,020 ranked candidates and take the first 1,000 with unique lift"
        ),
        "eligible_unique_canonical_snvs_all_splits": len(all_variants),
        "eligible_unique_canonical_snvs_gate_blocks": len(gate_variants),
        "target_variants": TARGET,
        "selected_candidates": len(rows),
        "selected_unique_1mb_blocks": len(block_counts),
        "largest_variants_per_block": max(block_counts.values()),
        "selected_by_chromosome": dict(sorted(chromosome_counts.items())),
        "sample_sha256": sample_sha256(rows),
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"candidates={len(rows)} target={TARGET} gate_eligible={len(gate_variants)} "
        f"blocks={len(block_counts)} sample_sha256={manifest['sample_sha256']}"
    )


if __name__ == "__main__":
    main()
