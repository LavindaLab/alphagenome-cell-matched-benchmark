#!/usr/bin/env python3
"""Publication-grade characterization of the frozen AlphaGenome confirmation result."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import sqlite3
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")
BUDGETS = (96, 384, 1000, 5000)
POLICIES = ("random", "pip", "avi", "magnitude", "concordance", "concordance_block_cap")
COVERAGES = (0.01, 0.05, 0.10, 0.25, 0.50)
BLOCK_SIZES = (1_000_000, 5_000_000)
SEED = 20260914


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, records) -> None:
    if not records:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_tsv_gz(path: Path, records) -> None:
    if not records:
        raise ValueError(f"refusing to write empty table: {path}")
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def load_frozen_module(project: Path):
    path = project / "tools/analyze_confirmation.py"
    spec = importlib.util.spec_from_file_location("frozen_confirmation_analysis", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def quantile(values, probability):
    return float(np.quantile(np.asarray(values, dtype=float), probability, method="linear"))


def interval(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return quantile(values, 0.025), quantile(values, 0.975)


def record_tuple(row):
    return row.chromosome, row.position, int(row.emvar)


def error_tuple(row, direction):
    return row.chromosome, row.position, int(direction != row.outcome_direction)


def bootstrap_rates(groups, block_size, replicates, seed):
    """Chromosome-stratified cluster bootstrap conditional on fixed groups."""
    keys = list(groups)
    all_blocks = sorted({(chrom, (pos - 1) // block_size) for values in groups.values() for chrom, pos, _ in values})
    block_index = {block: index for index, block in enumerate(all_blocks)}
    by_chromosome = defaultdict(list)
    for index, (chromosome, _) in enumerate(all_blocks):
        by_chromosome[chromosome].append(index)

    counts = np.zeros((len(keys), len(all_blocks)), dtype=np.float64)
    positives = np.zeros_like(counts)
    for group_index, key in enumerate(keys):
        for chromosome, position, outcome in groups[key]:
            index = block_index[(chromosome, (position - 1) // block_size)]
            counts[group_index, index] += 1
            positives[group_index, index] += outcome

    rng = np.random.default_rng(seed)
    weights = np.zeros((replicates, len(all_blocks)), dtype=np.int16)
    replicate_indices = np.arange(replicates)
    for indices in by_chromosome.values():
        choices = rng.choice(np.asarray(indices), size=(replicates, len(indices)), replace=True)
        for column in range(choices.shape[1]):
            np.add.at(weights, (replicate_indices, choices[:, column]), 1)
    numerators = positives @ weights.T
    denominators = counts @ weights.T
    with np.errstate(divide="ignore", invalid="ignore"):
        rates = numerators / denominators
    return {key: rates[index] for index, key in enumerate(keys)}


def average_precision_and_auroc(scores, outcomes):
    scores = np.asarray(scores, dtype=float)
    outcomes = np.asarray(outcomes, dtype=np.int8)
    positives = int(outcomes.sum())
    negatives = len(outcomes) - positives
    if not positives or not negatives:
        return math.nan, math.nan

    order = np.argsort(-scores, kind="mergesort")
    ordered_scores = scores[order]
    ordered_outcomes = outcomes[order]
    starts = np.r_[0, np.flatnonzero(ordered_scores[1:] != ordered_scores[:-1]) + 1]
    ends = np.r_[starts[1:], len(scores)]
    cumulative_positive = 0
    cumulative_total = 0
    average_precision = 0.0
    for start, end in zip(starts, ends, strict=True):
        group_positive = int(ordered_outcomes[start:end].sum())
        cumulative_positive += group_positive
        cumulative_total += end - start
        average_precision += (group_positive / positives) * (cumulative_positive / cumulative_total)

    ascending = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[ascending]
    sorted_outcomes = outcomes[ascending]
    starts = np.r_[0, np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]) + 1]
    ends = np.r_[starts[1:], len(scores)]
    positive_rank_sum = 0.0
    for start, end in zip(starts, ends, strict=True):
        mean_rank = ((start + 1) + end) / 2
        positive_rank_sum += mean_rank * int(sorted_outcomes[start:end].sum())
    auroc = (positive_rank_sum - positives * (positives + 1) / 2) / (positives * negatives)
    return float(average_precision), float(auroc)


def logistic_cluster(y, x, clusters):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    design = np.column_stack([np.ones(len(y)), x])
    beta = np.zeros(2, dtype=float)
    for _ in range(50):
        linear = np.clip(design @ beta, -30, 30)
        probability = 1 / (1 + np.exp(-linear))
        weights = np.clip(probability * (1 - probability), 1e-9, None)
        hessian = design.T @ (weights[:, None] * design)
        step = np.linalg.solve(hessian, design.T @ (y - probability))
        beta += step
        if float(np.max(np.abs(step))) < 1e-10:
            break
    linear = np.clip(design @ beta, -30, 30)
    probability = 1 / (1 + np.exp(-linear))
    weights = np.clip(probability * (1 - probability), 1e-9, None)
    bread = np.linalg.inv(design.T @ (weights[:, None] * design))
    scores = design * (y - probability)[:, None]
    unique, inverse = np.unique(np.asarray(clusters), return_inverse=True)
    cluster_scores = np.zeros((len(unique), 2), dtype=float)
    np.add.at(cluster_scores, inverse, scores)
    meat = cluster_scores.T @ cluster_scores
    correction = (len(unique) / (len(unique) - 1)) * ((len(y) - 1) / (len(y) - 2))
    covariance = correction * bread @ meat @ bread
    standard_error = math.sqrt(float(covariance[1, 1]))
    coefficient = float(beta[1])
    z_value = coefficient / standard_error
    p_value = math.erfc(abs(z_value) / math.sqrt(2))
    return coefficient, standard_error, p_value, len(unique)


def bh_adjust(records) -> None:
    eligible = [(index, float(row["p_value"])) for index, row in enumerate(records) if row["confirmation_role"] == "confirmatory"]
    ordered = sorted(eligible, key=lambda item: item[1])
    adjusted = [0.0] * len(ordered)
    running = 1.0
    for reverse_index in range(len(ordered) - 1, -1, -1):
        _, p_value = ordered[reverse_index]
        rank = reverse_index + 1
        running = min(running, p_value * len(ordered) / rank)
        adjusted[reverse_index] = running
    for (item, _), q_value in zip(ordered, adjusted, strict=True):
        records[item]["q_value_bh_20_tests"] = q_value
    for row in records:
        row.setdefault("q_value_bh_20_tests", "")


def verify_inputs(project, frozen):
    confirmation = project / "confirmation"
    results = confirmation / "results"
    authorization = read_json(confirmation / "OUTCOME_ACCESS_AUTHORIZATION.json")
    result_manifest = read_json(results / "result_manifest.json")
    inputs = {
        "metadata": confirmation / "variant_metadata.tsv.gz",
        "features": confirmation / "prediction_features.tsv.gz",
        "avi": confirmation / "avi_scores.tsv.gz",
        "regions": confirmation / "tss_distance.tsv.gz",
        "adastra": project / "evidence/adastra_exact_cell_confirmation.tsv.gz",
        "outcomes": confirmation / "mpra_outcomes_confirmation.tsv.gz",
        "outcome_receipt": confirmation / "outcome_extraction_receipt.json",
        "authorization": confirmation / "OUTCOME_ACCESS_AUTHORIZATION.json",
    }
    for name in ("metadata", "features", "avi", "regions", "adastra"):
        if sha256(inputs[name]) != authorization["analysis_input_sha256"][name]:
            raise ValueError(f"frozen input hash failed: {name}")
    if sha256(project / "tools/analyze_confirmation.py") != authorization["analysis_script_sha256"]:
        raise ValueError("frozen confirmatory analysis script changed")
    for name, expected in result_manifest["output_sha256"].items():
        if sha256(results / name) != expected:
            raise ValueError(f"frozen result hash failed: {name}")
    outcome_receipt = read_json(inputs["outcome_receipt"])
    if outcome_receipt["gate_outcomes_read"] is not False or outcome_receipt["scientific_rules_mutated"] is not False:
        raise ValueError("outcome firewall receipt failed")
    return inputs, authorization, result_manifest, outcome_receipt


def main() -> None:
    args = parse_args()
    project = args.project_dir.resolve()
    output = (args.output_dir or project / "manuscript/analysis").resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    frozen = load_frozen_module(project)
    inputs, authorization, result_manifest, outcome_receipt = verify_inputs(project, frozen)
    frozen_utility = read_tsv(project / "confirmation/results/budget_utility.tsv")
    utility_lookup = {(row["cell_context"], int(row["budget"]), row["policy"]): row for row in frozen_utility}

    descriptor, database_name = tempfile.mkstemp(prefix="alphagenome-manuscript-", suffix=".sqlite")
    os.close(descriptor)
    connection = sqlite3.connect(database_name)
    database_args = SimpleNamespace(
        metadata=inputs["metadata"], regions=inputs["regions"], avi=inputs["avi"],
        features=inputs["features"], outcomes=inputs["outcomes"],
    )
    import_counts = frozen.build_database(connection, database_args)

    budget_groups = {}
    risk_groups = {}
    rank_groups = {}
    rank_group_rows = []
    rank_metrics = []
    indicator_rows = []
    integrity_cells = {}
    try:
        for cell in CELLS:
            rows = frozen.load_cell(connection, cell)
            universe = [row for row in rows if row.fine_mapped]
            primary = [row for row in rows if row.enhancer_active and row.emvar]
            integrity_cells[cell] = {
                "measured_noncoding_variant_cells": len(rows),
                "eligible_fine_mapped": len(universe),
                "primary_enhancer_active_emvars": len(primary),
                "unique_variants": len({row.variant for row in rows}),
            }

            selections = {}
            for budget in BUDGETS:
                for policy in POLICIES:
                    selected = frozen.select_policy(universe, policy, budget, cell)
                    selections[(budget, policy)] = selected
                    budget_groups[(cell, budget, policy)] = [record_tuple(row) for row in selected]
                    expected = utility_lookup[(cell, budget, policy)]
                    observed_hits = sum(row.emvar for row in selected)
                    if len(selected) != int(expected["selected"]) or observed_hits != int(expected["emvar_hits"]):
                        raise ValueError(f"frozen selection replay failed: {cell} {budget} {policy}")

            for coverage in COVERAGES:
                requested = math.ceil(len(primary) * coverage)
                magnitude = sorted(primary, key=lambda row: (-row.magnitude, row.tie_hash))[:requested]
                concordance = sorted((row for row in primary if row.concordance_eligible), key=lambda row: (-row.concordance, row.tie_hash))[:requested]
                risk_groups[(cell, coverage, "magnitude")] = [error_tuple(row, row.magnitude_direction) for row in magnitude]
                risk_groups[(cell, coverage, "concordance")] = [error_tuple(row, row.concordance_direction) for row in concordance]

            requested = math.ceil(len(universe) * 0.10)
            magnitude_top = {row.variant for row in sorted(universe, key=lambda row: (-row.magnitude, row.tie_hash))[:requested]}
            avi_top = {row.variant for row in sorted(universe, key=lambda row: (-row.avi, row.tie_hash))[:requested]}
            grouped = defaultdict(list)
            for row in universe:
                in_magnitude = row.variant in magnitude_top
                in_avi = row.variant in avi_top
                group = "both_top10" if in_magnitude and in_avi else "magnitude_only_top10" if in_magnitude else "avi_only_top10" if in_avi else "neither_top10"
                grouped[group].append(row)
            for group in ("both_top10", "magnitude_only_top10", "avi_only_top10", "neither_top10"):
                values = grouped[group]
                rank_groups[(cell, group)] = [record_tuple(row) for row in values]
                rank_group_rows.append({
                    "cell_context": cell,
                    "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                    "group": group,
                    "variants": len(values),
                    "emvars": sum(row.emvar for row in values),
                    "emvar_prevalence": sum(row.emvar for row in values) / len(values) if values else math.nan,
                })

            for score_name, eligible, score_getter in (
                ("magnitude", universe, lambda row: row.magnitude),
                ("avi", universe, lambda row: row.avi),
                ("pip", [row for row in universe if row.pip is not None], lambda row: row.pip),
            ):
                ap, auc = average_precision_and_auroc([score_getter(row) for row in eligible], [row.emvar for row in eligible])
                rank_metrics.append({
                    "cell_context": cell,
                    "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                    "score": score_name,
                    "population": len(eligible),
                    "emvars": sum(row.emvar for row in eligible),
                    "prevalence": sum(row.emvar for row in eligible) / len(eligible),
                    "average_precision": ap,
                    "auroc": auc,
                })

            indicator_definitions = (
                ("dnase_h3k27ac_sign_disagreement", rows, np.asarray([row.sign_disagreement for row in rows], dtype=float)),
                ("absolute_quantile_difference_per_sd", rows, np.asarray([row.quantile_difference for row in rows], dtype=float)),
                ("dnase_exact_vs_other_direction_disagreement", rows, np.asarray([row.dnase_other_disagreement for row in rows], dtype=float)),
                ("h3k27ac_exact_vs_other_direction_disagreement", rows, np.asarray([row.histone_other_disagreement for row in rows], dtype=float)),
            )
            for name, population, predictor in indicator_definitions:
                if name.endswith("per_sd"):
                    predictor = (predictor - predictor.mean()) / predictor.std(ddof=0)
                y = np.asarray([row.emvar for row in population], dtype=float)
                clusters = [row.block_1mb for row in population]
                coefficient, standard_error, p_value, cluster_count = logistic_cluster(y, predictor, clusters)
                indicator_rows.append({
                    "cell_context": cell,
                    "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                    "indicator": name,
                    "population": len(population),
                    "clusters_1mb": cluster_count,
                    "odds_ratio": math.exp(coefficient),
                    "ci95_low": math.exp(coefficient - 1.959963984540054 * standard_error),
                    "ci95_high": math.exp(coefficient + 1.959963984540054 * standard_error),
                    "p_value": p_value,
                })
            discordant = grouped["avi_only_top10"] + grouped["magnitude_only_top10"]
            predictor = np.asarray([row.variant in magnitude_top for row in discordant], dtype=float)
            coefficient, standard_error, p_value, cluster_count = logistic_cluster(
                [row.emvar for row in discordant], predictor, [row.block_1mb for row in discordant]
            )
            indicator_rows.append({
                "cell_context": cell,
                "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                "indicator": "magnitude_only_vs_avi_only_top10",
                "population": len(discordant),
                "clusters_1mb": cluster_count,
                "odds_ratio": math.exp(coefficient),
                "ci95_low": math.exp(coefficient - 1.959963984540054 * standard_error),
                "ci95_high": math.exp(coefficient + 1.959963984540054 * standard_error),
                "p_value": p_value,
            })
            del rows, universe, primary
    finally:
        connection.close()
        Path(database_name).unlink(missing_ok=True)

    bh_adjust(indicator_rows)
    write_tsv(output / "ranking_metrics.tsv", rank_metrics)
    write_tsv(output / "uncertainty_indicator_models.tsv", indicator_rows)

    budget_contrasts = []
    budget_draw_records = []
    budget_rates_by_size = {}
    for block_size in BLOCK_SIZES:
        draws = bootstrap_rates(budget_groups, block_size, args.bootstrap_replicates, SEED + block_size)
        budget_rates_by_size[block_size] = draws
        for budget in BUDGETS:
            for comparator in ("random", "pip", "avi", "concordance", "concordance_block_cap"):
                cell_differences = []
                for cell in CELLS:
                    difference = draws[(cell, budget, "magnitude")] - draws[(cell, budget, comparator)]
                    cell_differences.append(difference)
                    low, high = interval(difference)
                    magnitude_row = utility_lookup[(cell, budget, "magnitude")]
                    comparator_row = utility_lookup[(cell, budget, comparator)]
                    point = float(magnitude_row["precision"]) - float(comparator_row["precision"])
                    budget_contrasts.append({
                        "cell_context": cell,
                        "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                        "budget": budget,
                        "comparator": comparator,
                        "block_size": block_size,
                        "magnitude_precision": magnitude_row["precision"],
                        "comparator_precision": comparator_row["precision"],
                        "precision_difference": point,
                        "ci95_low": low,
                        "ci95_high": high,
                        "magnitude_hits": magnitude_row["emvar_hits"],
                        "comparator_hits": comparator_row["emvar_hits"],
                        "additional_hits": int(magnitude_row["emvar_hits"]) - int(comparator_row["emvar_hits"]),
                    })
                confirmation_macro = np.mean([cell_differences[CELLS.index(cell)] for cell in CONFIRMATION_CELLS], axis=0)
                point_values = []
                magnitude_hits = comparator_hits = 0
                for cell in CONFIRMATION_CELLS:
                    magnitude_row = utility_lookup[(cell, budget, "magnitude")]
                    comparator_row = utility_lookup[(cell, budget, comparator)]
                    point_values.append(float(magnitude_row["precision"]) - float(comparator_row["precision"]))
                    magnitude_hits += int(magnitude_row["emvar_hits"])
                    comparator_hits += int(comparator_row["emvar_hits"])
                low, high = interval(confirmation_macro)
                budget_contrasts.append({
                    "cell_context": "confirmation_macro",
                    "confirmation_role": "confirmatory_macro",
                    "budget": budget,
                    "comparator": comparator,
                    "block_size": block_size,
                    "magnitude_precision": np.mean([float(utility_lookup[(cell, budget, "magnitude")]["precision"]) for cell in CONFIRMATION_CELLS]),
                    "comparator_precision": np.mean([float(utility_lookup[(cell, budget, comparator)]["precision"]) for cell in CONFIRMATION_CELLS]),
                    "precision_difference": np.mean(point_values),
                    "ci95_low": low,
                    "ci95_high": high,
                    "magnitude_hits": magnitude_hits,
                    "comparator_hits": comparator_hits,
                    "additional_hits": magnitude_hits - comparator_hits,
                })
            central_differences = {
                cell: draws[(cell, budget, "magnitude")] - draws[(cell, budget, "avi")] for cell in CELLS
            }
            for replicate in range(args.bootstrap_replicates):
                budget_draw_records.append({
                    "block_size": block_size,
                    "replicate": replicate + 1,
                    "budget": budget,
                    "confirmation_macro_difference": np.mean([central_differences[cell][replicate] for cell in CONFIRMATION_CELLS]),
                    **{f"{cell}_difference": central_differences[cell][replicate] for cell in CELLS},
                })
    write_tsv(output / "budget_precision_contrasts.tsv", budget_contrasts)
    write_tsv_gz(output / "budget_precision_bootstrap_draws.tsv.gz", budget_draw_records)

    risk_intervals = []
    risk_contrasts = []
    for block_size in BLOCK_SIZES:
        draws = bootstrap_rates(risk_groups, block_size, args.bootstrap_replicates, SEED + block_size + 31)
        for coverage in COVERAGES:
            differences = []
            for cell in CELLS:
                for policy in ("magnitude", "concordance"):
                    values = draws[(cell, coverage, policy)]
                    low, high = interval(values)
                    risk_intervals.append({
                        "cell_context": cell,
                        "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                        "coverage": coverage,
                        "policy": policy,
                        "block_size": block_size,
                        "sign_error": float(np.mean([value[2] for value in risk_groups[(cell, coverage, policy)]])),
                        "ci95_low": low,
                        "ci95_high": high,
                        "selected": len(risk_groups[(cell, coverage, policy)]),
                    })
                difference = draws[(cell, coverage, "concordance")] - draws[(cell, coverage, "magnitude")]
                differences.append(difference)
                low, high = interval(difference)
                risk_contrasts.append({
                    "cell_context": cell,
                    "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                    "coverage": coverage,
                    "block_size": block_size,
                    "concordance_minus_magnitude": float(
                        np.mean([value[2] for value in risk_groups[(cell, coverage, "concordance")]])
                        - np.mean([value[2] for value in risk_groups[(cell, coverage, "magnitude")]])
                    ),
                    "ci95_low": low,
                    "ci95_high": high,
                })
            macro = np.mean([differences[CELLS.index(cell)] for cell in CONFIRMATION_CELLS], axis=0)
            low, high = interval(macro)
            risk_contrasts.append({
                "cell_context": "confirmation_macro",
                "confirmation_role": "confirmatory_macro",
                "coverage": coverage,
                "block_size": block_size,
                "concordance_minus_magnitude": np.mean([
                    np.mean([value[2] for value in risk_groups[(cell, coverage, "concordance")]])
                    - np.mean([value[2] for value in risk_groups[(cell, coverage, "magnitude")]])
                    for cell in CONFIRMATION_CELLS
                ]),
                "ci95_low": low,
                "ci95_high": high,
            })
    write_tsv(output / "selective_risk_intervals.tsv", risk_intervals)
    write_tsv(output / "selective_risk_contrasts.tsv", risk_contrasts)

    rank_intervals = []
    rank_contrasts = []
    rank_rates_by_size = {}
    for block_size in BLOCK_SIZES:
        draws = bootstrap_rates(rank_groups, block_size, args.bootstrap_replicates, SEED + block_size + 73)
        rank_rates_by_size[block_size] = draws
        for cell in CELLS:
            for group in ("both_top10", "magnitude_only_top10", "avi_only_top10", "neither_top10"):
                values = draws[(cell, group)]
                low, high = interval(values)
                exact = next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == group)
                rank_intervals.append({**exact, "block_size": block_size, "ci95_low": low, "ci95_high": high})
            magnitude_only = draws[(cell, "magnitude_only_top10")]
            avi_only = draws[(cell, "avi_only_top10")]
            difference = magnitude_only - avi_only
            ratio = magnitude_only / avi_only
            difference_low, difference_high = interval(difference)
            ratio_low, ratio_high = interval(ratio)
            mag_exact = next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "magnitude_only_top10")
            avi_exact = next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "avi_only_top10")
            rank_contrasts.append({
                "cell_context": cell,
                "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                "block_size": block_size,
                "magnitude_only_prevalence": mag_exact["emvar_prevalence"],
                "avi_only_prevalence": avi_exact["emvar_prevalence"],
                "prevalence_difference": float(mag_exact["emvar_prevalence"]) - float(avi_exact["emvar_prevalence"]),
                "difference_ci95_low": difference_low,
                "difference_ci95_high": difference_high,
                "prevalence_ratio": float(mag_exact["emvar_prevalence"]) / float(avi_exact["emvar_prevalence"]),
                "ratio_ci95_low": ratio_low,
                "ratio_ci95_high": ratio_high,
            })
        macro_difference = np.mean([
            draws[(cell, "magnitude_only_top10")] - draws[(cell, "avi_only_top10")] for cell in CONFIRMATION_CELLS
        ], axis=0)
        macro_ratio = np.mean([
            draws[(cell, "magnitude_only_top10")] / draws[(cell, "avi_only_top10")] for cell in CONFIRMATION_CELLS
        ], axis=0)
        low, high = interval(macro_difference)
        ratio_low, ratio_high = interval(macro_ratio)
        exact_mag = np.mean([
            float(next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "magnitude_only_top10")["emvar_prevalence"])
            for cell in CONFIRMATION_CELLS
        ])
        exact_avi = np.mean([
            float(next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "avi_only_top10")["emvar_prevalence"])
            for cell in CONFIRMATION_CELLS
        ])
        rank_contrasts.append({
            "cell_context": "confirmation_macro",
            "confirmation_role": "confirmatory_macro",
            "block_size": block_size,
            "magnitude_only_prevalence": exact_mag,
            "avi_only_prevalence": exact_avi,
            "prevalence_difference": exact_mag - exact_avi,
            "difference_ci95_low": low,
            "difference_ci95_high": high,
            "prevalence_ratio": np.mean([
                float(next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "magnitude_only_top10")["emvar_prevalence"])
                / float(next(row for row in rank_group_rows if row["cell_context"] == cell and row["group"] == "avi_only_top10")["emvar_prevalence"])
                for cell in CONFIRMATION_CELLS
            ]),
            "ratio_ci95_low": ratio_low,
            "ratio_ci95_high": ratio_high,
        })
    write_tsv(output / "avi_context_rank_groups.tsv", rank_intervals)
    write_tsv(output / "avi_context_rank_contrasts.tsv", rank_contrasts)

    central_1mb = [row for row in budget_contrasts if row["block_size"] == 1_000_000 and row["comparator"] == "avi"]
    central_5mb = [row for row in budget_contrasts if row["block_size"] == 5_000_000 and row["comparator"] == "avi"]
    rank_central = [row for row in rank_contrasts if row["cell_context"] == "confirmation_macro"]
    summary = {
        "analysis_role": "post-confirmation manuscript characterization; not a new confirmatory analysis",
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_seed": SEED,
        "fixed_budgets": list(BUDGETS),
        "fixed_coverages": list(COVERAGES),
        "confirmation_cells": list(CONFIRMATION_CELLS),
        "a549_role": "descriptive_count_rule_downgrade",
        "magnitude_minus_avi": {
            str(budget): {
                "one_mb": next(row for row in central_1mb if row["budget"] == budget and row["cell_context"] == "confirmation_macro"),
                "five_mb": next(row for row in central_5mb if row["budget"] == budget and row["cell_context"] == "confirmation_macro"),
            } for budget in BUDGETS
        },
        "rank_disagreement_top10_confirmation_macro": {
            "one_mb": next(row for row in rank_central if row["block_size"] == 1_000_000),
            "five_mb": next(row for row in rank_central if row["block_size"] == 5_000_000),
        },
        "frozen_primary": read_json(project / "confirmation/results/confirmation_results.json"),
        "scientific_rules_mutated": False,
        "new_threshold_optimized_on_outcomes": False,
        "gate_outcomes_used": False,
    }
    (output / "manuscript_analysis_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    reference = read_json(project / "confirmation/reference_validation_hg38_summary.json")
    atlas_inspection = read_json(project / "confirmation/atlas/inspection_receipt.json")
    retrieval = read_json(project / "confirmation/atlas/retrieval_summary.json")
    result_inspection = read_json(project / "confirmation/results/inspection_receipt.json")
    integrity = {
        "all_checks_passed": True,
        "protocol_sha256": sha256(project / "protocol/PROTOCOL_v0.1.md"),
        "protocol_sha256_authorized": authorization["protocol_sha256"],
        "evidence_ledger_sha256": sha256(project / "protocol/EVIDENCE_LEDGER.md"),
        "evidence_ledger_sha256_authorized": authorization["evidence_ledger_sha256"],
        "frozen_analysis_sha256": sha256(project / "tools/analyze_confirmation.py"),
        "frozen_analysis_sha256_authorized": authorization["analysis_script_sha256"],
        "result_checks": f"{result_inspection['checks_passed']}/{result_inspection['checks_total']}",
        "atlas_checks": f"{atlas_inspection['checks_passed']}/{atlas_inspection['checks_total']}",
        "atlas_complete": retrieval["complete"],
        "atlas_incomplete": retrieval["incomplete"],
        "atlas_request_attempts": retrieval["request_attempts"],
        "atlas_retries": retrieval["retries"],
        "exact_ontology_terms": retrieval["ontology_terms"],
        "reference_validation_pre_exclusion": reference,
        "outcome_extraction_after_freeze": outcome_receipt["confirmation_outcomes_extracted_after_freeze"],
        "gate_outcomes_read": outcome_receipt["gate_outcomes_read"],
        "scientific_rules_mutated": False,
        "database_import_counts": import_counts,
        "cell_counts_and_uniqueness": integrity_cells,
        "duplicate_variant_cell_rows": sum(value["measured_noncoding_variant_cells"] - value["unique_variants"] for value in integrity_cells.values()),
    }
    expected_protocol = integrity["protocol_sha256"] == integrity["protocol_sha256_authorized"]
    expected_evidence = integrity["evidence_ledger_sha256"] == integrity["evidence_ledger_sha256_authorized"]
    if not expected_protocol or not expected_evidence or integrity["duplicate_variant_cell_rows"] != 0:
        raise ValueError("integrity audit failed")
    (output / "integrity_audit.json").write_text(json.dumps(integrity, indent=2, sort_keys=True) + "\n")

    output_hashes = {path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()}
    manifest = {
        "analysis_script_sha256": sha256(Path(__file__)),
        "frozen_result_manifest_sha256": sha256(project / "confirmation/results/result_manifest.json"),
        "output_sha256": output_hashes,
        "bootstrap_replicates": args.bootstrap_replicates,
        "scientific_rules_mutated": False,
        "gate_outcomes_used": False,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    core = summary["magnitude_minus_avi"]
    print(
        "passed "
        f"budget384={core['384']['one_mb']['precision_difference']:.6f} "
        f"ci=[{core['384']['one_mb']['ci95_low']:.6f},{core['384']['one_mb']['ci95_high']:.6f}] "
        f"budget1000={core['1000']['one_mb']['precision_difference']:.6f} "
        f"ci=[{core['1000']['one_mb']['ci95_low']:.6f},{core['1000']['one_mb']['ci95_high']:.6f}]"
    )


if __name__ == "__main__":
    main()
