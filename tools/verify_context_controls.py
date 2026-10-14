#!/usr/bin/env python3
"""Independently recompute and verify the post hoc context-control headlines."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")
BUDGETS = (384, 1000)
SCORES = ("dnase_only", "h3k27ac_only", "combined_magnitude", "avi")
BLOCK_SIZES = (1_000_000, 5_000_000)
REPLICATES = 2000
BASE_SEED = 20261008
TOLERANCE = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def truth(value: str) -> int:
    return int(value == "true")


def batch_insert(database, statement, values, batch_size=10000):
    batch = []
    count = 0
    for value in values:
        batch.append(value)
        if len(batch) == batch_size:
            database.executemany(statement, batch)
            count += len(batch)
            batch.clear()
    if batch:
        database.executemany(statement, batch)
        count += len(batch)
    database.commit()
    return count


def close(left, right, tolerance=TOLERANCE):
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=tolerance)


def ci(values):
    values = np.asarray(values, dtype=float)
    return (
        float(np.quantile(values, 0.025, method="linear")),
        float(np.quantile(values, 0.975, method="linear")),
    )


def build_database(database, project):
    confirmation = project / "confirmation"
    database.executescript("""
        PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE; PRAGMA cache_size=-500000;
        CREATE TABLE metadata(variant TEXT PRIMARY KEY, fine INTEGER);
        CREATE TABLE regions(variant TEXT PRIMARY KEY);
        CREATE TABLE avi(variant TEXT PRIMARY KEY, score REAL);
        CREATE TABLE outcomes(variant TEXT, source TEXT, cell TEXT, emvar INTEGER, PRIMARY KEY(variant,cell));
        CREATE TABLE features(
          variant TEXT, cell TEXT, magnitude REAL, concordance REAL, winner TEXT, tie_hash TEXT,
          PRIMARY KEY(variant,cell));
    """)
    counts = {}
    counts["metadata"] = batch_insert(
        database,
        "INSERT INTO metadata VALUES (?,?)",
        ((row["variant_id_hg38"], truth(row["fine_mapped_candidate"])) for row in rows(confirmation / "variant_metadata.tsv.gz")),
    )
    counts["regions"] = batch_insert(
        database,
        "INSERT INTO regions VALUES (?)",
        ((row["variant_id_hg38"],) for row in rows(confirmation / "tss_distance.tsv.gz")),
    )
    counts["avi"] = batch_insert(
        database,
        "INSERT INTO avi VALUES (?,?)",
        ((row["variant_id_hg38"], float(row["quantile_score"])) for row in rows(confirmation / "avi_scores.tsv.gz")),
    )
    counts["outcomes"] = batch_insert(
        database,
        "INSERT INTO outcomes VALUES (?,?,?,?)",
        (
            (row["variant_id_hg38"], row["source_variant_id_hg19"], row["cell_context"], truth(row["emvar"]))
            for row in rows(confirmation / "mpra_outcomes_confirmation.tsv.gz")
            if row["measured"] == "true" and row["cell_context"] in CONFIRMATION_CELLS
        ),
    )
    identity = {
        "rows_checked": 0,
        "max_failures": 0,
        "min_failures": 0,
        "winner_failures": 0,
        "ties_at_tolerance": 0,
        "row_level_sample": None,
    }

    def feature_values():
        for row in rows(confirmation / "prediction_features.tsv.gz"):
            magnitude = float(row["magnitude_M"])
            concordance = float(row["concordance_C"])
            winner = row["magnitude_winning_modality"]
            if winner == "DNASE":
                dnase, h3k27ac = magnitude, concordance
            elif winner == "H3K27ac":
                dnase, h3k27ac = concordance, magnitude
            else:
                dnase = h3k27ac = math.nan
            identity["rows_checked"] += 1
            identity["ties_at_tolerance"] += int(abs(magnitude - concordance) <= TOLERANCE)
            identity["max_failures"] += int(not close(max(dnase, h3k27ac), magnitude))
            identity["min_failures"] += int(not close(min(dnase, h3k27ac), concordance))
            assigned = dnase if winner == "DNASE" else h3k27ac
            identity["winner_failures"] += int(winner not in ("DNASE", "H3K27ac") or not close(assigned, magnitude))
            if identity["row_level_sample"] is None:
                identity["row_level_sample"] = {
                    "variant_id_hg38": row["variant_id_hg38"],
                    "cell_context": row["cell_context"],
                    "magnitude_M": magnitude,
                    "concordance_C": concordance,
                    "winner": winner,
                    "reconstructed_abs_dnase_quantile": dnase,
                    "reconstructed_abs_h3k27ac_quantile": h3k27ac,
                }
            yield (
                row["variant_id_hg38"],
                row["cell_context"],
                magnitude,
                concordance,
                winner,
                row["tie_break_sha256"],
            )

    counts["features"] = batch_insert(database, "INSERT INTO features VALUES (?,?,?,?,?,?)", feature_values())
    database.executescript("CREATE INDEX outcomes_cell ON outcomes(cell); CREATE INDEX features_variant ON features(variant);")
    return counts, identity


def candidates_for_cell(database, outcome_cell):
    query = """
      SELECT o.variant,o.source,o.emvar,a.score,f.cell,f.magnitude,f.concordance,f.winner,f.tie_hash
      FROM outcomes o
      JOIN metadata m ON m.variant=o.variant AND m.fine=1
      JOIN regions r ON r.variant=o.variant
      JOIN avi a ON a.variant=o.variant
      JOIN features f ON f.variant=o.variant
      WHERE o.cell=? ORDER BY o.variant,f.cell
    """
    candidates = []
    current = None
    base = None
    by_cell = {}
    tie_hashes = set()

    def finish():
        if base is None:
            return
        if set(by_cell) != set(CELLS) or len(tie_hashes) != 1:
            raise ValueError(f"incomplete or inconsistent prediction contexts for {base[0]}")
        magnitude, concordance, winner = by_cell[outcome_cell]
        dnase, h3k27ac = (magnitude, concordance) if winner == "DNASE" else (concordance, magnitude)
        chromosome, position, _, _ = base[1].split(":")
        candidates.append({
            "variant": base[0],
            "chromosome": chromosome,
            "position": int(position),
            "emvar": int(base[2]),
            "avi": float(base[3]),
            "tie_hash": next(iter(tie_hashes)),
            "combined_magnitude": magnitude,
            "dnase_only": dnase,
            "h3k27ac_only": h3k27ac,
            "magnitude_by_cell": {cell: values[0] for cell, values in by_cell.items()},
        })

    for row in database.execute(query, (outcome_cell,)):
        if current is not None and row[0] != current:
            finish()
            by_cell = {}
            tie_hashes = set()
        current = row[0]
        base = row[:4]
        by_cell[row[4]] = (float(row[5]), float(row[6]), row[7])
        tie_hashes.add(row[8])
    finish()
    return candidates


def select(candidates, budget, score_getter):
    return sorted(candidates, key=lambda row: (-score_getter(row), row["tie_hash"]))[:budget]


def selected_values(selected):
    return [(row["chromosome"], row["position"], row["emvar"]) for row in selected]


def independent_bootstrap(groups, block_size, replicates, seed):
    """Independent implementation of the fixed-selection stratified block bootstrap."""
    keys = list(groups)
    blocks = sorted({
        (chromosome, (position - 1) // block_size)
        for values in groups.values()
        for chromosome, position, _ in values
    })
    index = {block: number for number, block in enumerate(blocks)}
    chromosome_indices = defaultdict(list)
    for number, (chromosome, _) in enumerate(blocks):
        chromosome_indices[chromosome].append(number)
    counts = np.zeros((len(keys), len(blocks)), dtype=np.float64)
    positives = np.zeros_like(counts)
    for group_number, key in enumerate(keys):
        for chromosome, position, outcome in groups[key]:
            block_number = index[(chromosome, (position - 1) // block_size)]
            counts[group_number, block_number] += 1
            positives[group_number, block_number] += outcome
    rng = np.random.default_rng(seed)
    weights = np.zeros((replicates, len(blocks)), dtype=np.int16)
    replicate_numbers = np.arange(replicates)
    for indices in chromosome_indices.values():
        sampled = rng.choice(np.asarray(indices), size=(replicates, len(indices)), replace=True)
        for column in range(sampled.shape[1]):
            np.add.at(weights, (replicate_numbers, sampled[:, column]), 1)
    numerators = positives @ weights.T
    denominators = counts @ weights.T
    with np.errstate(divide="ignore", invalid="ignore"):
        rates = numerators / denominators
    return {key: rates[number] for number, key in enumerate(keys)}


def main() -> None:
    args = parse_args()
    project = args.project_dir.resolve()
    output = (args.output_dir or project / "manuscript/analysis/context-controls").resolve()
    required = (
        "single_track_estimates.tsv",
        "single_track_contrasts.tsv",
        "wrong_cell_grid.tsv",
        "wrong_cell_contrasts.tsv",
        "context_control_summary.json",
    )
    for name in required:
        if not (output / name).is_file():
            raise ValueError(f"missing result file: {name}")

    observed_single_estimates = list(rows(output / "single_track_estimates.tsv"))
    observed_single_contrasts = list(rows(output / "single_track_contrasts.tsv"))
    observed_wrong_grid = list(rows(output / "wrong_cell_grid.tsv"))
    observed_wrong_contrasts = list(rows(output / "wrong_cell_contrasts.tsv"))
    summary = json.loads((output / "context_control_summary.json").read_text())

    checks = []

    def check(name, condition, observed=None, expected=None):
        checks.append({
            "check": name,
            "passed": bool(condition),
            "observed": observed,
            "expected": expected,
        })

    descriptor, database_name = tempfile.mkstemp(prefix="verify-context-controls-", suffix=".sqlite")
    os.close(descriptor)
    database = sqlite3.connect(database_name)
    groups = {}
    point = {}
    eligibility = {}
    try:
        import_counts, identity = build_database(database, project)
        check("reconstruction_max_identity_all_rows", identity["max_failures"] == 0, identity["max_failures"], 0)
        check("reconstruction_min_identity_all_rows", identity["min_failures"] == 0, identity["min_failures"], 0)
        check("reconstruction_winner_identity_all_rows", identity["winner_failures"] == 0, identity["winner_failures"], 0)
        for outcome_cell in CONFIRMATION_CELLS:
            candidates = candidates_for_cell(database, outcome_cell)
            eligibility[outcome_cell] = len(candidates)
            for budget in BUDGETS:
                for score in SCORES:
                    selected_rows = select(candidates, budget, lambda row, name=score: row[name])
                    key = ("single", outcome_cell, budget, score)
                    groups[key] = selected_values(selected_rows)
                    point[key] = (len(selected_rows), sum(row["emvar"] for row in selected_rows))
                for prediction_cell in CELLS:
                    selected_rows = select(
                        candidates,
                        budget,
                        lambda row, cell=prediction_cell: row["magnitude_by_cell"][cell],
                    )
                    key = ("context", outcome_cell, budget, prediction_cell)
                    groups[key] = selected_values(selected_rows)
                    point[key] = (len(selected_rows), sum(row["emvar"] for row in selected_rows))
    finally:
        database.close()
        Path(database_name).unlink(missing_ok=True)

    check("single_track_table_row_count", len(observed_single_estimates) == 80, len(observed_single_estimates), 80)
    check("single_track_contrast_row_count", len(observed_single_contrasts) == 100, len(observed_single_contrasts), 100)
    check("wrong_cell_grid_row_count", len(observed_wrong_grid) == 80, len(observed_wrong_grid), 80)
    check("wrong_cell_contrast_row_count", len(observed_wrong_contrasts) == 20, len(observed_wrong_contrasts), 20)
    check("all_rows_explicitly_post_hoc", all(
        row["post_hoc"] == "true"
        for table in (observed_single_estimates, observed_single_contrasts, observed_wrong_grid, observed_wrong_contrasts)
        for row in table
    ))
    check("summary_explicitly_post_hoc", summary["analysis_status"] == "post_hoc", summary["analysis_status"], "post_hoc")
    check("four_by_five_grid_complete", {
        (row["outcome_cell"], row["prediction_cell"], int(row["budget"]), int(row["block_size_bp"]))
        for row in observed_wrong_grid
    } == {
        (outcome, prediction, budget, block)
        for outcome in CONFIRMATION_CELLS
        for prediction in CELLS
        for budget in BUDGETS
        for block in BLOCK_SIZES
    })

    point_failures = []
    for row in observed_single_estimates:
        budget = int(row["budget"])
        score = row["score"]
        if row["cell_context"] in CONFIRMATION_CELLS:
            selected_count, hits = point[("single", row["cell_context"], budget, score)]
            expected_precision = hits / selected_count
        else:
            values = [point[("single", cell, budget, score)] for cell in CONFIRMATION_CELLS]
            selected_count = sum(value[0] for value in values)
            hits = sum(value[1] for value in values)
            expected_precision = float(np.mean([value[1] / value[0] for value in values]))
        if not (
            int(row["selected"]) == selected_count
            and int(row["emvar_hits"]) == hits
            and close(row["precision"], expected_precision)
        ):
            point_failures.append([row["cell_context"], budget, score])
    check("all_single_track_points_recomputed", not point_failures, point_failures, [])

    grid_failures = []
    for row in observed_wrong_grid:
        key = ("context", row["outcome_cell"], int(row["budget"]), row["prediction_cell"])
        selected_count, hits = point[key]
        if not (
            int(row["selected"]) == selected_count
            and int(row["emvar_hits"]) == hits
            and close(row["precision"], hits / selected_count)
        ):
            grid_failures.append(list(key[1:]))
    check("all_wrong_cell_grid_points_recomputed", not grid_failures, grid_failures, [])

    authoritative = {
        (row["cell_context"], int(row["budget"]), row["policy"]): row
        for row in rows(project / "confirmation/results/budget_utility.tsv")
        if row["cell_context"] in CONFIRMATION_CELLS and int(row["budget"]) in BUDGETS
    }
    replay_failures = []
    for cell in CONFIRMATION_CELLS:
        for budget in BUDGETS:
            for score, policy in (("combined_magnitude", "magnitude"), ("avi", "avi")):
                selected_count, hits = point[("single", cell, budget, score)]
                expected = authoritative[(cell, budget, policy)]
                if selected_count != int(expected["selected"]) or hits != int(expected["emvar_hits"]):
                    replay_failures.append([cell, budget, score])
    check("combined_and_avi_authoritative_replay", not replay_failures, replay_failures, [])

    headline_failures = []
    comparisons = (
        ("dnase_only", "combined_magnitude"),
        ("h3k27ac_only", "combined_magnitude"),
        ("dnase_only", "avi"),
        ("h3k27ac_only", "avi"),
        ("combined_magnitude", "avi"),
    )
    for block_size in BLOCK_SIZES:
        draws = independent_bootstrap(groups, block_size, REPLICATES, BASE_SEED + block_size)
        for budget in BUDGETS:
            for focal, reference in comparisons:
                macro = np.mean([
                    draws[("single", cell, budget, focal)] - draws[("single", cell, budget, reference)]
                    for cell in CONFIRMATION_CELLS
                ], axis=0)
                expected_low, expected_high = ci(macro)
                observed = next(
                    row for row in observed_single_contrasts
                    if row["cell_context"] == "confirmation_macro"
                    and int(row["budget"]) == budget
                    and int(row["block_size_bp"]) == block_size
                    and row["score"] == focal
                    and row["reference_score"] == reference
                )
                focal_points = [point[("single", cell, budget, focal)] for cell in CONFIRMATION_CELLS]
                reference_points = [point[("single", cell, budget, reference)] for cell in CONFIRMATION_CELLS]
                expected_difference = float(np.mean([hits / selected for selected, hits in focal_points]) - np.mean([hits / selected for selected, hits in reference_points]))
                if not (
                    close(observed["precision_difference"], expected_difference)
                    and close(observed["ci95_low"], expected_low)
                    and close(observed["ci95_high"], expected_high)
                ):
                    headline_failures.append([block_size, budget, focal, reference])
            cell_differences = []
            point_differences = []
            for outcome_cell in CONFIRMATION_CELLS:
                exact = ("context", outcome_cell, budget, outcome_cell)
                wrong = [
                    ("context", outcome_cell, budget, prediction)
                    for prediction in CELLS
                    if prediction != outcome_cell
                ]
                cell_differences.append(draws[exact] - np.mean([draws[key] for key in wrong], axis=0))
                exact_selected, exact_hits = point[exact]
                point_differences.append(exact_hits / exact_selected - np.mean([point[key][1] / point[key][0] for key in wrong]))
            macro = np.mean(cell_differences, axis=0)
            expected_low, expected_high = ci(macro)
            expected_difference = float(np.mean(point_differences))
            observed = next(
                row for row in observed_wrong_contrasts
                if row["cell_context"] == "confirmation_macro"
                and int(row["budget"]) == budget
                and int(row["block_size_bp"]) == block_size
            )
            if not (
                close(observed["exact_minus_all_wrong_mean"], expected_difference)
                and close(observed["ci95_low"], expected_low)
                and close(observed["ci95_high"], expected_high)
            ):
                headline_failures.append([block_size, budget, "exact", "all_wrong_mean"])
    check("all_macro_headlines_and_intervals_independently_recomputed", not headline_failures, headline_failures, [])

    numeric_fields = (
        (observed_single_estimates, ("precision", "ci95_low", "ci95_high")),
        (observed_single_contrasts, ("precision_difference", "ci95_low", "ci95_high")),
        (observed_wrong_grid, ("precision", "ci95_low", "ci95_high")),
        (observed_wrong_contrasts, ("exact_minus_all_wrong_mean", "ci95_low", "ci95_high")),
    )
    check("all_reported_estimates_finite", all(
        math.isfinite(float(row[field]))
        for table, fields in numeric_fields
        for row in table
        for field in fields
    ))

    failed = [item["check"] for item in checks if not item["passed"]]
    receipt = {
        "verification_role": "independent recomputation from released derived inputs; post hoc controls",
        "passed": not failed,
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "failed_checks": failed,
        "identity_validation": identity,
        "database_import_counts": import_counts,
        "eligible_variants_by_outcome_cell": eligibility,
        "bootstrap_replicates": REPLICATES,
        "bootstrap_seed_by_block_size": {str(size): BASE_SEED + size for size in BLOCK_SIZES},
        "checks": checks,
    }
    (output / "verification_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    if failed:
        raise ValueError(f"context-control verification failed: {failed}")
    print(f"passed checks={len(checks)}/{len(checks)} identities={identity['rows_checked']} headline_failures=0")


if __name__ == "__main__":
    main()
