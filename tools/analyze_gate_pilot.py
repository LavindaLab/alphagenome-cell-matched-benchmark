#!/usr/bin/env python3
"""Apply the frozen primary and experimental-utility policies to gate outcomes."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import random
import statistics
from datetime import datetime, timezone
from pathlib import Path


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
COVERAGES = (0.01, 0.05, 0.10, 0.25, 0.50)
BUDGETS = (96, 384, 1000, 5000)
BOOTSTRAPS = 2000
SEED = 20260913


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def truth(raw: str) -> bool:
    return raw == "true"


def sign(value: float) -> int:
    return (value > 0) - (value < 0)


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def evaluate_primary(rows: list[dict[str, object]], coverage: float):
    requested = math.ceil(len(rows) * coverage)
    magnitude = sorted(rows, key=lambda row: (-float(row["magnitude_M"]), row["tie_break_sha256"]))[:requested]
    concordant = sorted(
        (row for row in rows if row["concordance_eligible"]),
        key=lambda row: (-float(row["concordance_C"]), row["tie_break_sha256"]),
    )[:requested]
    def error(selected, direction_field):
        return statistics.fmean(
            int(int(row[direction_field]) != int(row["outcome_direction"])) for row in selected
        ) if selected else math.nan
    magnitude_error = error(magnitude, "magnitude_predicted_direction")
    concordance_error = error(concordant, "concordance_predicted_direction")
    return {
        "population": len(rows),
        "requested": requested,
        "magnitude_selected": len(magnitude),
        "concordance_selected": len(concordant),
        "magnitude_error": magnitude_error,
        "concordance_error": concordance_error,
        "difference": concordance_error - magnitude_error,
    }


def bootstrap(primary_by_cell, block_size: int, seed: int, continuation_cells):
    rng = random.Random(seed)
    union_blocks = collections.defaultdict(set)
    for rows in primary_by_cell.values():
        for row in rows:
            chromosome, position, _, _ = str(row["source_variant_id_hg19"]).split(":")
            union_blocks[chromosome].add((int(position) - 1) // block_size)
    draws = []
    for replicate in range(BOOTSTRAPS):
        multiplicity = collections.Counter()
        for chromosome, blocks in sorted(union_blocks.items()):
            ordered = sorted(blocks)
            for _ in ordered:
                multiplicity[(chromosome, rng.choice(ordered))] += 1
        differences = []
        cell_differences = {}
        for cell in CELLS:
            sampled = []
            for row in primary_by_cell[cell]:
                chromosome, position, _, _ = str(row["source_variant_id_hg19"]).split(":")
                copies = multiplicity[(chromosome, (int(position) - 1) // block_size)]
                sampled.extend([row] * copies)
            estimate = evaluate_primary(sampled, 0.10)
            if math.isfinite(estimate["difference"]):
                differences.append(estimate["difference"])
                cell_differences[cell] = estimate["difference"]
        continuation_differences = [
            cell_differences[cell] for cell in continuation_cells
            if cell in cell_differences
        ]
        draws.append({
            "replicate": replicate + 1,
            "block_size": block_size,
            "macro_difference": statistics.fmean(differences),
            "continuation_macro_difference": statistics.fmean(
                continuation_differences
            ),
            **{
                f"{cell}_difference": cell_differences.get(cell, "")
                for cell in CELLS
            },
        })
    return draws


def random_stratified(rows, budget):
    strata = collections.defaultdict(list)
    for row in rows:
        chromosome = str(row["source_variant_id_hg19"]).split(":", 1)[0]
        strata[(chromosome, bool(row["pip_available"]))].append(row)
    total = len(rows)
    quotas = {key: budget * len(values) / total for key, values in strata.items()}
    allocation = {key: math.floor(value) for key, value in quotas.items()}
    remaining = budget - sum(allocation.values())
    order = sorted(strata, key=lambda key: (-(quotas[key] - allocation[key]), str(key)))
    for key in order[:remaining]:
        allocation[key] += 1
    selected = []
    for key, values in strata.items():
        selected.extend(sorted(values, key=lambda row: row["random_sha256"])[:allocation[key]])
    return selected


def select_budget(rows, policy: str, budget: int):
    requested = min(budget, len(rows))
    if policy == "random":
        return random_stratified(rows, requested)
    if policy == "pip":
        eligible = [row for row in rows if row["max_pip"] is not None]
        return sorted(eligible, key=lambda row: (-float(row["max_pip"]), row["tie_break_sha256"]))[:requested]
    if policy == "avi":
        eligible = [row for row in rows if row["avi_score"] is not None]
        return sorted(eligible, key=lambda row: (-float(row["avi_score"]), row["tie_break_sha256"]))[:requested]
    if policy == "magnitude":
        return sorted(rows, key=lambda row: (-float(row["magnitude_M"]), row["tie_break_sha256"]))[:requested]
    concordant = sorted(
        (row for row in rows if row["concordance_eligible"]),
        key=lambda row: (-float(row["concordance_C"]), row["tie_break_sha256"]),
    )
    if policy == "concordance":
        return concordant[:requested]
    if policy == "concordance_block_cap":
        selected = []
        block_counts = collections.Counter()
        for row in concordant:
            block = row["block_id_1mb_hg19"]
            if block_counts[block] >= 2:
                continue
            selected.append(row)
            block_counts[block] += 1
            if len(selected) == requested:
                break
        return selected
    raise ValueError(policy)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("outcomes", type=Path)
    parser.add_argument("variant_metadata", type=Path)
    parser.add_argument("noncoding", type=Path)
    parser.add_argument("prediction_features", type=Path)
    parser.add_argument("avi_scores", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outcomes = read_tsv(args.outcomes)
    metadata = {row["variant_id_hg38"]: row for row in read_tsv(args.variant_metadata)}
    noncoding = {row["variant_id_hg38"]: truth(row["noncoding_eligible"]) for row in read_tsv(args.noncoding)}
    features = {(row["variant_id_hg38"], row["cell_context"]): row for row in read_tsv(args.prediction_features)}
    manifest = {
        row["variant_id"]: row
        for row in read_tsv(args.outcomes.parent / "variants_hg38.tsv")
    }
    avi = {
        row["variant_id_hg38"]: float(row["quantile_score"])
        for row in (json.loads(line) for line in args.avi_scores.read_text().splitlines())
    }
    analysis_rows = []
    for outcome in outcomes:
        key = (outcome["variant_id_hg38"], outcome["cell_context"])
        feature = features.get(key)
        if feature is None or not noncoding[outcome["variant_id_hg38"]] or not truth(outcome["measured"]):
            continue
        variant_meta = metadata[outcome["variant_id_hg38"]]
        manifest_row = manifest[outcome["variant_id_hg38"]]
        source = outcome["source_variant_id_hg19"]
        analysis_rows.append({
            **outcome,
            "active": truth(outcome["active"]),
            "enhancer_active": truth(outcome["enhancer_active"]),
            "emvar": truth(outcome["emvar"]),
            "outcome_direction": sign(float(outcome["log2skew_alt_minus_ref"])),
            "magnitude_M": float(feature["magnitude_M"]),
            "magnitude_predicted_direction": int(feature["magnitude_predicted_direction"]),
            "concordance_eligible": truth(feature["concordance_eligible"]),
            "concordance_C": float(feature["concordance_C"]),
            "concordance_predicted_direction": int(feature["concordance_predicted_direction"]),
            "sign_disagreement": truth(feature["dnase_h3k27ac_sign_disagreement"]),
            "tie_break_sha256": feature["tie_break_sha256"],
            "random_sha256": hashlib.sha256(f"AGATLAS-UTILITY-RANDOM-2026-V01\0{source}\0{outcome['cell_context']}".encode()).hexdigest(),
            "fine_mapped_candidate": truth(variant_meta["fine_mapped_candidate"]),
            "pip_available": truth(variant_meta["pip_available"]),
            "max_pip": None if variant_meta["max_pip"] == "" else float(variant_meta["max_pip"]),
            "avi_score": avi.get(outcome["variant_id_hg38"]),
            "sampling_weight": float(manifest_row["stratum_population"])
            / float(manifest_row["stratum_target"]),
        })

    primary_by_cell = {
        cell: [row for row in analysis_rows if row["cell_context"] == cell and row["enhancer_active"] and row["emvar"]]
        for cell in CELLS
    }
    primary_rows = []
    risk_rows = []
    for cell in CELLS:
        for coverage in COVERAGES:
            estimate = evaluate_primary(primary_by_cell[cell], coverage)
            risk_rows.append({"cell_context": cell, "coverage": coverage, **estimate})
            if coverage == 0.10:
                primary_rows.append({"cell_context": cell, **estimate})
    macro_difference = statistics.fmean(row["difference"] for row in primary_rows)
    continuation_cells = [
        row["cell_context"] for row in primary_rows if row["population"] >= 200
    ]
    downgraded_cells = [
        row["cell_context"] for row in primary_rows if row["population"] < 200
    ]
    continuation_macro_difference = statistics.fmean(
        row["difference"] for row in primary_rows
        if row["cell_context"] in continuation_cells
    )
    bootstrap_1mb = bootstrap(
        primary_by_cell, 1_000_000, SEED, continuation_cells
    )
    bootstrap_5mb = bootstrap(
        primary_by_cell, 5_000_000, SEED + 5, continuation_cells
    )
    bootstrap_rows = bootstrap_1mb + bootstrap_5mb
    values_1mb = [row["macro_difference"] for row in bootstrap_1mb]
    values_5mb = [row["macro_difference"] for row in bootstrap_5mb]
    continuation_values_1mb = [
        row["continuation_macro_difference"] for row in bootstrap_1mb
    ]
    continuation_values_5mb = [
        row["continuation_macro_difference"] for row in bootstrap_5mb
    ]

    precision_rows = []
    for cell in CELLS:
        cell_rows = [row for row in analysis_rows if row["cell_context"] == cell]
        prevalence = statistics.fmean(int(row["emvar"]) for row in cell_rows)
        for coverage in COVERAGES:
            requested = math.ceil(len(cell_rows) * coverage)
            for policy in ("magnitude", "concordance"):
                selected = select_budget(cell_rows, policy, requested)
                precision = statistics.fmean(int(row["emvar"]) for row in selected) if selected else math.nan
                precision_rows.append({
                    "cell_context": cell, "coverage": coverage, "policy": policy,
                    "eligible_population": len(cell_rows), "requested": requested,
                    "selected": len(selected), "emvar_hits": sum(row["emvar"] for row in selected),
                    "precision": precision, "prevalence": prevalence,
                    "precision_lift": precision / prevalence if prevalence else math.nan,
                    "missed_emvar_fraction": (sum(row["emvar"] for row in cell_rows) - sum(row["emvar"] for row in selected)) / sum(row["emvar"] for row in cell_rows) if any(row["emvar"] for row in cell_rows) else math.nan,
                })

    utility_rows = []
    policies = ("random", "pip", "avi", "magnitude", "concordance", "concordance_block_cap")
    for cell in CELLS:
        universe = [row for row in analysis_rows if row["cell_context"] == cell and row["fine_mapped_candidate"]]
        for budget in BUDGETS:
            for policy in policies:
                selected = select_budget(universe, policy, budget)
                utility_rows.append({
                    "cell_context": cell, "budget": budget, "policy": policy,
                    "eligible_population": len(universe), "selected": len(selected),
                    "emvar_hits": sum(row["emvar"] for row in selected),
                    "ipw_emvar_hits": sum(float(row["sampling_weight"]) for row in selected if row["emvar"]),
                    "independent_blocks_with_hit": len({row["block_id_1mb_hg19"] for row in selected if row["emvar"]}),
                    "enhancer_active_selected": sum(row["enhancer_active"] for row in selected),
                    "avi_missing_selected_universe": sum(row["avi_score"] is None for row in universe),
                })

    discordance_rows = []
    for cell in CELLS:
        for disagrees in (False, True):
            subset = [row for row in analysis_rows if row["cell_context"] == cell and row["sign_disagreement"] == disagrees]
            emvars = [row for row in subset if row["emvar"] and row["enhancer_active"]]
            discordance_rows.append({
                "cell_context": cell, "dnase_h3k27ac_sign_disagreement": str(disagrees).lower(),
                "measured_noncoding": len(subset), "emvars": sum(row["emvar"] for row in subset),
                "emvar_prevalence": statistics.fmean(int(row["emvar"]) for row in subset) if subset else math.nan,
                "enhancer_active_emvars": len(emvars),
                "magnitude_sign_error_among_enhancer_active_emvars": statistics.fmean(int(row["magnitude_predicted_direction"] != row["outcome_direction"]) for row in emvars) if emvars else math.nan,
            })

    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "frozen 10,000-SNV gate-block pilot only",
        "primary_estimand": "five-cell macro mean sign-error(concordant)-sign-error(magnitude) at 10% coverage among noncoding enhancer-active emVars",
        "primary_macro_difference": macro_difference,
        "primary_1mb_bootstrap_95_ci": [percentile(values_1mb, 0.025), percentile(values_1mb, 0.975)],
        "primary_5mb_bootstrap_95_ci": [percentile(values_5mb, 0.025), percentile(values_5mb, 0.975)],
        "continuation_primary_cells": continuation_cells,
        "downgraded_cells_under_frozen_200_count_rule": downgraded_cells,
        "continuation_macro_difference": continuation_macro_difference,
        "continuation_1mb_bootstrap_95_ci": [
            percentile(continuation_values_1mb, 0.025),
            percentile(continuation_values_1mb, 0.975),
        ],
        "continuation_5mb_bootstrap_95_ci": [
            percentile(continuation_values_5mb, 0.025),
            percentile(continuation_values_5mb, 0.975),
        ],
        "primary_counts_by_cell": {row["cell_context"]: row["population"] for row in primary_rows},
        "cells_meeting_200_emvar_continuation_count": continuation_cells,
        "cells_below_200_emvar_continuation_count": downgraded_cells,
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": SEED,
        "confirmation_split_accessed": False,
        "scientific_rules_mutated": False,
    }
    write_tsv(args.output_dir / "primary_by_cell.tsv", primary_rows)
    write_tsv(args.output_dir / "risk_coverage.tsv", risk_rows)
    write_tsv(args.output_dir / "precision_lift.tsv", precision_rows)
    write_tsv(args.output_dir / "budget_utility.tsv", utility_rows)
    write_tsv(args.output_dir / "discordance.tsv", discordance_rows)
    write_tsv(args.output_dir / "bootstrap_draws.tsv", bootstrap_rows)
    (args.output_dir / "pilot_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"primary_difference={macro_difference:.6f} "
        f"ci95=[{result['primary_1mb_bootstrap_95_ci'][0]:.6f},{result['primary_1mb_bootstrap_95_ci'][1]:.6f}] "
        f"counts={result['primary_counts_by_cell']}"
    )


if __name__ == "__main__":
    main()
