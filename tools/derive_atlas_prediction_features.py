#!/usr/bin/env python3
"""Exercise frozen prediction-only policy plumbing without MPRA labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def sign(value: float) -> int:
    return (value > 0) - (value < 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("atlas_scores_jsonl", type=Path)
    parser.add_argument("features_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    args = parser.parse_args()

    scores = [
        json.loads(line) for line in args.atlas_scores_jsonl.read_text().splitlines()
    ]
    grouped: dict[tuple[str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    by_variant_modality: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for row in scores:
        key = (str(row["variant_id_hg38"]), str(row["cell_context"]))
        grouped[key][str(row["modality"])] = row
        by_variant_modality[
            (str(row["variant_id_hg38"]), str(row["modality"]))
        ][str(row["cell_context"])] = float(row["raw_score"])

    rows = []
    for (variant, cell), modalities in sorted(grouped.items()):
        if set(modalities) != {"DNASE", "H3K27ac"}:
            raise ValueError(f"incomplete modality pair for {variant} {cell}")
        dnase = modalities["DNASE"]
        histone = modalities["H3K27ac"]
        raw_dnase = float(dnase["raw_score"])
        raw_histone = float(histone["raw_score"])
        q_dnase = float(dnase["quantile_score"])
        q_histone = float(histone["quantile_score"])
        magnitude_modality = (
            "DNASE" if abs(q_dnase) >= abs(q_histone) else "H3K27ac"
        )
        magnitude_raw = raw_dnase if magnitude_modality == "DNASE" else raw_histone
        sign_dnase = sign(raw_dnase)
        sign_histone = sign(raw_histone)
        concordant = sign_dnase != 0 and sign_dnase == sign_histone
        other_disagreement = {}
        for modality, exact_raw in (
            ("DNASE", raw_dnase),
            ("H3K27ac", raw_histone),
        ):
            other_values = [
                value
                for other_cell, value in by_variant_modality[(variant, modality)].items()
                if other_cell != cell
            ]
            other_sign = sign(statistics.median(other_values))
            other_disagreement[modality] = sign(exact_raw) != other_sign
        rows.append(
            {
                "variant_id_hg38": variant,
                "cell_context": cell,
                "magnitude_M": max(abs(q_dnase), abs(q_histone)),
                "magnitude_winning_modality": magnitude_modality,
                "magnitude_predicted_direction": sign(magnitude_raw),
                "concordance_eligible": str(concordant).lower(),
                "concordance_C": min(abs(q_dnase), abs(q_histone)),
                "concordance_predicted_direction": sign_dnase if concordant else 0,
                "dnase_h3k27ac_sign_disagreement": str(
                    sign_dnase != sign_histone
                ).lower(),
                "absolute_quantile_difference": abs(q_dnase - q_histone),
                "dnase_exact_vs_other_median_direction_disagreement": str(
                    other_disagreement["DNASE"]
                ).lower(),
                "h3k27ac_exact_vs_other_median_direction_disagreement": str(
                    other_disagreement["H3K27ac"]
                ).lower(),
                "tie_break_sha256": hashlib.sha256(variant.encode()).hexdigest(),
            }
        )
    if not all(
        math.isfinite(float(row[field]))
        for row in rows
        for field in ("magnitude_M", "concordance_C", "absolute_quantile_difference")
    ):
        raise ValueError("non-finite derived feature")

    with args.features_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    cell_counts = Counter(str(row["cell_context"]) for row in rows)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "variant_cell_rows": len(rows),
        "unique_variants": len({row["variant_id_hg38"] for row in rows}),
        "rows_by_cell": dict(sorted(cell_counts.items())),
        "definitions_exercised": [
            "frozen magnitude M and winning-modality direction with DNase ties",
            "frozen concordance eligibility, C and common-sign direction",
            "named prediction-only discordance indicators",
            "frozen SHA-256 tie-break field",
        ],
        "coverage_threshold_applied": False,
        "ranking_selection_applied": False,
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"prediction_features=pass rows={len(rows)} variants={summary['unique_variants']}")


if __name__ == "__main__":
    main()
