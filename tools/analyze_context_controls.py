#!/usr/bin/env python3
"""Run the bounded post hoc single-track and cell-context controls."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import sqlite3
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")
BUDGETS = (384, 1000)
SCORES = ("dnase_only", "h3k27ac_only", "combined_magnitude", "avi")
BLOCK_SIZES = (1_000_000, 5_000_000)
BOOTSTRAPS = 2000
SEED = 20261008
TOLERANCE = 1e-12
POST_HOC_DATE = "2026-10-08"
ORIGINAL_COMMIT = "d4ac41b345b69bb6f1ec5324699f9f7e3c96bcc7"


@dataclass(slots=True)
class Candidate:
    variant: str
    source: str
    chromosome: str
    position: int
    block_1mb: str
    emvar: bool
    avi: float
    tie_hash: str
    magnitude: float
    dnase_only: float
    h3k27ac_only: float
    magnitude_by_cell: dict[str, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAPS)
    return parser.parse_args()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read_json(path: Path):
    return json.loads(path.read_text())


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, records) -> None:
    records = list(records)
    if not records:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_tsv_gz_deterministic(path: Path, records) -> None:
    records = list(records)
    if not records:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with _gzip_text_writer(compressed) as handle:
                writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(records)


class _gzip_text_writer:
    def __init__(self, compressed):
        self.handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")

    def __enter__(self):
        return self.handle

    def __exit__(self, exc_type, exc_value, traceback):
        self.handle.flush()
        self.handle.detach()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def interval(values) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        raise ValueError("bootstrap interval contains no finite values")
    return (
        float(np.quantile(array, 0.025, method="linear")),
        float(np.quantile(array, 0.975, method="linear")),
    )


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def verify_frozen_inputs(project: Path, frozen) -> dict[str, Path]:
    confirmation = project / "confirmation"
    results = confirmation / "results"
    paths = {
        "metadata": confirmation / "variant_metadata.tsv.gz",
        "features": confirmation / "prediction_features.tsv.gz",
        "avi": confirmation / "avi_scores.tsv.gz",
        "regions": confirmation / "tss_distance.tsv.gz",
        "outcomes": confirmation / "mpra_outcomes_confirmation.tsv.gz",
        "authorization": confirmation / "OUTCOME_ACCESS_AUTHORIZATION.json",
        "outcome_receipt": confirmation / "outcome_extraction_receipt.json",
        "result_manifest": results / "result_manifest.json",
        "utility": results / "budget_utility.tsv",
    }
    authorization = read_json(paths["authorization"])
    for name in ("metadata", "features", "avi", "regions"):
        if sha256(paths[name]) != authorization["analysis_input_sha256"][name]:
            raise ValueError(f"frozen input hash failed: {name}")
    if sha256(project / "tools/analyze_confirmation.py") != authorization["analysis_script_sha256"]:
        raise ValueError("frozen confirmatory analysis script hash failed")
    outcome_receipt = read_json(paths["outcome_receipt"])
    if sha256(paths["outcomes"]) != outcome_receipt["output_sha256"]:
        raise ValueError("released outcome hash failed")
    result_manifest = read_json(paths["result_manifest"])
    for name, expected in result_manifest["output_sha256"].items():
        if sha256(results / name) != expected:
            raise ValueError(f"authoritative result hash failed: {name}")
    if outcome_receipt["gate_outcomes_read"] is not False:
        raise ValueError("outcome receipt does not preserve the gate/confirmation firewall")
    return paths


def create_database(database: sqlite3.Connection, paths, frozen) -> tuple[dict, dict]:
    database.executescript("""
        PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE; PRAGMA cache_size=-500000;
        CREATE TABLE metadata(variant TEXT PRIMARY KEY, fine INTEGER);
        CREATE TABLE regions(variant TEXT PRIMARY KEY);
        CREATE TABLE avi(variant TEXT PRIMARY KEY, score REAL);
        CREATE TABLE outcomes(
          variant TEXT, source TEXT, cell TEXT, block_1mb TEXT, emvar INTEGER,
          PRIMARY KEY(variant, cell));
        CREATE TABLE features(
          variant TEXT, cell TEXT, magnitude REAL, concordance REAL, winner TEXT,
          sign_disagreement INTEGER, absolute_difference REAL, tie_hash TEXT,
          PRIMARY KEY(variant, cell));
    """)
    counts = {}
    counts["metadata"] = frozen.batch_insert(
        database,
        "INSERT INTO metadata VALUES (?,?)",
        ((row["variant_id_hg38"], frozen.truth(row["fine_mapped_candidate"])) for row in frozen.rows(paths["metadata"])),
    )
    counts["regions"] = frozen.batch_insert(
        database,
        "INSERT INTO regions VALUES (?)",
        ((row["variant_id_hg38"],) for row in frozen.rows(paths["regions"])),
    )
    counts["avi"] = frozen.batch_insert(
        database,
        "INSERT INTO avi VALUES (?,?)",
        ((row["variant_id_hg38"], float(row["quantile_score"])) for row in frozen.rows(paths["avi"])),
    )
    counts["outcomes"] = frozen.batch_insert(
        database,
        "INSERT INTO outcomes VALUES (?,?,?,?,?)",
        (
            (
                row["variant_id_hg38"],
                row["source_variant_id_hg19"],
                row["cell_context"],
                row["block_id_1mb_hg19"],
                frozen.truth(row["emvar"]),
            )
            for row in frozen.rows(paths["outcomes"])
            if row["measured"] == "true" and row["cell_context"] in CONFIRMATION_CELLS
        ),
    )

    validation = {
        "rows_checked": 0,
        "invalid_rows": 0,
        "ties_at_tolerance": 0,
        "max_identity_failures": 0,
        "min_identity_failures": 0,
        "winner_identity_failures": 0,
        "tolerance": TOLERANCE,
        "row_level_sample": None,
    }

    def feature_values():
        for row in frozen.rows(paths["features"]):
            magnitude = float(row["magnitude_M"])
            concordance = float(row["concordance_C"])
            winner = row["magnitude_winning_modality"]
            disagreement = frozen.truth(row["dnase_h3k27ac_sign_disagreement"])
            observed_difference = float(row["absolute_quantile_difference"])
            if winner == "DNASE":
                dnase, h3k27ac = magnitude, concordance
            elif winner == "H3K27ac":
                dnase, h3k27ac = concordance, magnitude
            else:
                validation["winner_identity_failures"] += 1
                dnase = h3k27ac = math.nan
            validation["rows_checked"] += 1
            if abs(magnitude - concordance) <= TOLERANCE:
                validation["ties_at_tolerance"] += 1
            max_ok = math.isclose(max(dnase, h3k27ac), magnitude, rel_tol=0.0, abs_tol=TOLERANCE)
            min_ok = math.isclose(min(dnase, h3k27ac), concordance, rel_tol=0.0, abs_tol=TOLERANCE)
            winner_value = dnase if winner == "DNASE" else h3k27ac
            winner_ok = math.isclose(winner_value, magnitude, rel_tol=0.0, abs_tol=TOLERANCE)
            validation["max_identity_failures"] += int(not max_ok)
            validation["min_identity_failures"] += int(not min_ok)
            validation["winner_identity_failures"] += int(not winner_ok)
            if not (max_ok and min_ok and winner_ok):
                validation["invalid_rows"] += 1
            if validation["row_level_sample"] is None:
                validation["row_level_sample"] = {
                    "variant_id_hg38": row["variant_id_hg38"],
                    "cell_context": row["cell_context"],
                    "magnitude_M": magnitude,
                    "concordance_C": concordance,
                    "magnitude_winning_modality": winner,
                    "reconstructed_abs_dnase_quantile": dnase,
                    "reconstructed_abs_h3k27ac_quantile": h3k27ac,
                    "reconstructed_max": max(dnase, h3k27ac),
                    "reconstructed_min": min(dnase, h3k27ac),
                }
            yield (
                row["variant_id_hg38"],
                row["cell_context"],
                magnitude,
                concordance,
                winner,
                disagreement,
                observed_difference,
                row["tie_break_sha256"],
            )

    counts["features"] = frozen.batch_insert(
        database,
        "INSERT INTO features VALUES (?,?,?,?,?,?,?,?)",
        feature_values(),
    )
    database.executescript("""
        CREATE INDEX outcomes_cell ON outcomes(cell);
        CREATE INDEX features_variant ON features(variant);
    """)
    if validation["invalid_rows"]:
        raise ValueError(f"single-track reconstruction identity failed for {validation['invalid_rows']} rows")
    return counts, validation


def load_candidates(database: sqlite3.Connection, outcome_cell: str) -> list[Candidate]:
    query = """
      SELECT o.variant,o.source,o.block_1mb,o.emvar,a.score,
             f.cell,f.magnitude,f.concordance,f.winner,f.tie_hash
      FROM outcomes o
      JOIN metadata m ON m.variant=o.variant AND m.fine=1
      JOIN regions r ON r.variant=o.variant
      JOIN avi a ON a.variant=o.variant
      JOIN features f ON f.variant=o.variant
      WHERE o.cell=?
      ORDER BY o.variant,f.cell
    """
    candidates = []
    current_variant = None
    base = None
    feature_values = {}
    tie_hashes = set()

    def finish_group():
        if base is None:
            return
        if set(feature_values) != set(CELLS):
            raise ValueError(f"candidate {base[0]} lacks complete five-context prediction coverage")
        if len(tie_hashes) != 1:
            raise ValueError(f"candidate {base[0]} has context-dependent tie hashes")
        exact_magnitude, exact_concordance, exact_winner = feature_values[outcome_cell]
        if exact_winner == "DNASE":
            dnase, h3k27ac = exact_magnitude, exact_concordance
        elif exact_winner == "H3K27ac":
            dnase, h3k27ac = exact_concordance, exact_magnitude
        else:
            raise ValueError(f"candidate {base[0]} has invalid winner")
        chromosome, position, _, _ = base[1].split(":")
        candidates.append(Candidate(
            variant=base[0],
            source=base[1],
            chromosome=chromosome,
            position=int(position),
            block_1mb=base[2],
            emvar=bool(base[3]),
            avi=float(base[4]),
            tie_hash=next(iter(tie_hashes)),
            magnitude=exact_magnitude,
            dnase_only=dnase,
            h3k27ac_only=h3k27ac,
            magnitude_by_cell={cell: values[0] for cell, values in feature_values.items()},
        ))

    for row in database.execute(query, (outcome_cell,)):
        variant = row[0]
        if current_variant is not None and variant != current_variant:
            finish_group()
            feature_values = {}
            tie_hashes = set()
        current_variant = variant
        base = row[:5]
        feature_values[row[5]] = (float(row[6]), float(row[7]), row[8])
        tie_hashes.add(row[9])
    finish_group()
    return candidates


def select_novel_score(candidates, budget, getter):
    return sorted(candidates, key=lambda row: (-getter(row), row.tie_hash))[: min(budget, len(candidates))]


def group_values(selected):
    return [(row.chromosome, row.position, int(row.emvar)) for row in selected]


def point_record(selected):
    hits = sum(row.emvar for row in selected)
    return {
        "selected": len(selected),
        "emvar_hits": hits,
        "precision": hits / len(selected) if selected else math.nan,
    }


def main() -> None:
    args = parse_args()
    if args.bootstrap_replicates != BOOTSTRAPS:
        raise ValueError(f"this settled analysis uses the repository-standard {BOOTSTRAPS} bootstrap replicates")
    project = args.project_dir.resolve()
    output = (args.output_dir or project / "manuscript/analysis/context-controls").resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    frozen = load_module(project / "tools/analyze_confirmation.py", "context_control_frozen_analysis")
    manuscript = load_module(project / "tools/analyze_manuscript.py", "context_control_manuscript_analysis")
    paths = verify_frozen_inputs(project, frozen)
    authoritative = read_tsv(paths["utility"])
    authoritative_lookup = {
        (row["cell_context"], int(row["budget"]), row["policy"]): row
        for row in authoritative
        if row["cell_context"] in CONFIRMATION_CELLS and int(row["budget"]) in BUDGETS
    }

    descriptor, database_name = tempfile.mkstemp(prefix="alphagenome-context-controls-", suffix=".sqlite")
    os.close(descriptor)
    database = sqlite3.connect(database_name)
    groups = {}
    points = {}
    eligibility = {}
    authoritative_checks = []
    try:
        import_counts, reconstruction = create_database(database, paths, frozen)
        for outcome_cell in CONFIRMATION_CELLS:
            candidates = load_candidates(database, outcome_cell)
            eligibility[outcome_cell] = {
                "eligible_variants": len(candidates),
                "prediction_contexts": list(CELLS),
                "identical_candidate_set_all_prediction_contexts": True,
            }
            expected_population = int(authoritative_lookup[(outcome_cell, BUDGETS[0], "magnitude")]["eligible_population"])
            if len(candidates) != expected_population:
                raise ValueError(
                    f"eligible population differs from authoritative result: {outcome_cell} "
                    f"observed={len(candidates)} expected={expected_population}"
                )
            for budget in BUDGETS:
                selected_by_score = {
                    "combined_magnitude": frozen.select_policy(candidates, "magnitude", budget, outcome_cell),
                    "avi": frozen.select_policy(candidates, "avi", budget, outcome_cell),
                    "dnase_only": select_novel_score(candidates, budget, lambda row: row.dnase_only),
                    "h3k27ac_only": select_novel_score(candidates, budget, lambda row: row.h3k27ac_only),
                }
                for score in SCORES:
                    selected = selected_by_score[score]
                    key = ("single", outcome_cell, budget, score)
                    groups[key] = group_values(selected)
                    points[key] = point_record(selected)
                for score, policy in (("combined_magnitude", "magnitude"), ("avi", "avi")):
                    observed = points[("single", outcome_cell, budget, score)]
                    expected = authoritative_lookup[(outcome_cell, budget, policy)]
                    passed = (
                        observed["selected"] == int(expected["selected"])
                        and observed["emvar_hits"] == int(expected["emvar_hits"])
                        and math.isclose(observed["precision"], float(expected["precision"]), rel_tol=0.0, abs_tol=1e-15)
                    )
                    authoritative_checks.append({
                        "cell_context": outcome_cell,
                        "budget": budget,
                        "score": score,
                        "observed_selected": observed["selected"],
                        "expected_selected": int(expected["selected"]),
                        "observed_hits": observed["emvar_hits"],
                        "expected_hits": int(expected["emvar_hits"]),
                        "observed_precision": observed["precision"],
                        "expected_precision": float(expected["precision"]),
                        "passed": passed,
                    })
                    if not passed:
                        raise ValueError(f"authoritative selection replay failed: {outcome_cell} {budget} {score}")
                for prediction_cell in CELLS:
                    selected = select_novel_score(
                        candidates,
                        budget,
                        lambda row, cell=prediction_cell: row.magnitude_by_cell[cell],
                    )
                    key = ("context", outcome_cell, budget, prediction_cell)
                    groups[key] = group_values(selected)
                    points[key] = point_record(selected)
    finally:
        database.close()
        Path(database_name).unlink(missing_ok=True)

    single_estimates = []
    single_contrasts = []
    wrong_grid = []
    wrong_contrasts = []
    bootstrap_records = []
    single_comparisons = (
        ("dnase_only", "combined_magnitude"),
        ("h3k27ac_only", "combined_magnitude"),
        ("dnase_only", "avi"),
        ("h3k27ac_only", "avi"),
        ("combined_magnitude", "avi"),
    )
    headline = {str(budget): {} for budget in BUDGETS}

    for block_size in BLOCK_SIZES:
        bootstrap_seed = SEED + block_size
        draws = manuscript.bootstrap_rates(groups, block_size, args.bootstrap_replicates, bootstrap_seed)
        for budget in BUDGETS:
            for score in SCORES:
                cell_draws = []
                for cell in CONFIRMATION_CELLS:
                    key = ("single", cell, budget, score)
                    low, high = interval(draws[key])
                    point = points[key]
                    single_estimates.append({
                        "post_hoc": "true",
                        "cell_context": cell,
                        "summary_level": "cell",
                        "budget": budget,
                        "score": score,
                        "block_size_bp": block_size,
                        **point,
                        "ci95_low": low,
                        "ci95_high": high,
                        "bootstrap_replicates": args.bootstrap_replicates,
                        "bootstrap_seed": bootstrap_seed,
                    })
                    cell_draws.append(draws[key])
                macro_draw = np.mean(cell_draws, axis=0)
                low, high = interval(macro_draw)
                macro_points = [points[("single", cell, budget, score)] for cell in CONFIRMATION_CELLS]
                single_estimates.append({
                    "post_hoc": "true",
                    "cell_context": "confirmation_macro",
                    "summary_level": "unweighted_four_cell_macro",
                    "budget": budget,
                    "score": score,
                    "block_size_bp": block_size,
                    "selected": sum(row["selected"] for row in macro_points),
                    "emvar_hits": sum(row["emvar_hits"] for row in macro_points),
                    "precision": float(np.mean([row["precision"] for row in macro_points])),
                    "ci95_low": low,
                    "ci95_high": high,
                    "bootstrap_replicates": args.bootstrap_replicates,
                    "bootstrap_seed": bootstrap_seed,
                })

            for focal, reference in single_comparisons:
                difference_draws = []
                for cell in CONFIRMATION_CELLS:
                    focal_key = ("single", cell, budget, focal)
                    reference_key = ("single", cell, budget, reference)
                    difference = draws[focal_key] - draws[reference_key]
                    difference_draws.append(difference)
                    low, high = interval(difference)
                    focal_point = points[focal_key]
                    reference_point = points[reference_key]
                    single_contrasts.append({
                        "post_hoc": "true",
                        "cell_context": cell,
                        "summary_level": "cell",
                        "budget": budget,
                        "score": focal,
                        "reference_score": reference,
                        "block_size_bp": block_size,
                        "score_precision": focal_point["precision"],
                        "reference_precision": reference_point["precision"],
                        "precision_difference": focal_point["precision"] - reference_point["precision"],
                        "ci95_low": low,
                        "ci95_high": high,
                        "score_hits": focal_point["emvar_hits"],
                        "reference_hits": reference_point["emvar_hits"],
                        "additional_hits": focal_point["emvar_hits"] - reference_point["emvar_hits"],
                        "bootstrap_replicates": args.bootstrap_replicates,
                        "bootstrap_seed": bootstrap_seed,
                    })
                macro_difference = np.mean(difference_draws, axis=0)
                low, high = interval(macro_difference)
                focal_points = [points[("single", cell, budget, focal)] for cell in CONFIRMATION_CELLS]
                reference_points = [points[("single", cell, budget, reference)] for cell in CONFIRMATION_CELLS]
                single_contrasts.append({
                    "post_hoc": "true",
                    "cell_context": "confirmation_macro",
                    "summary_level": "unweighted_four_cell_macro",
                    "budget": budget,
                    "score": focal,
                    "reference_score": reference,
                    "block_size_bp": block_size,
                    "score_precision": float(np.mean([row["precision"] for row in focal_points])),
                    "reference_precision": float(np.mean([row["precision"] for row in reference_points])),
                    "precision_difference": float(np.mean([row["precision"] for row in focal_points]) - np.mean([row["precision"] for row in reference_points])),
                    "ci95_low": low,
                    "ci95_high": high,
                    "score_hits": sum(row["emvar_hits"] for row in focal_points),
                    "reference_hits": sum(row["emvar_hits"] for row in reference_points),
                    "additional_hits": sum(row["emvar_hits"] for row in focal_points) - sum(row["emvar_hits"] for row in reference_points),
                    "bootstrap_replicates": args.bootstrap_replicates,
                    "bootstrap_seed": bootstrap_seed,
                })

            wrong_cell_differences = []
            exact_cell_draws = []
            wrong_mean_draws = []
            for outcome_cell in CONFIRMATION_CELLS:
                for prediction_cell in CELLS:
                    key = ("context", outcome_cell, budget, prediction_cell)
                    low, high = interval(draws[key])
                    relationship = "exact" if outcome_cell == prediction_cell else "wrong"
                    wrong_grid.append({
                        "post_hoc": "true",
                        "outcome_cell": outcome_cell,
                        "prediction_cell": prediction_cell,
                        "relationship": relationship,
                        "budget": budget,
                        "block_size_bp": block_size,
                        **points[key],
                        "ci95_low": low,
                        "ci95_high": high,
                        "bootstrap_replicates": args.bootstrap_replicates,
                        "bootstrap_seed": bootstrap_seed,
                    })
                exact_key = ("context", outcome_cell, budget, outcome_cell)
                wrong_keys = [
                    ("context", outcome_cell, budget, prediction_cell)
                    for prediction_cell in CELLS
                    if prediction_cell != outcome_cell
                ]
                wrong_draw = np.mean([draws[key] for key in wrong_keys], axis=0)
                difference = draws[exact_key] - wrong_draw
                wrong_cell_differences.append(difference)
                exact_cell_draws.append(draws[exact_key])
                wrong_mean_draws.append(wrong_draw)
                low, high = interval(difference)
                exact_point = points[exact_key]
                wrong_points = [points[key] for key in wrong_keys]
                wrong_contrasts.append({
                    "post_hoc": "true",
                    "cell_context": outcome_cell,
                    "summary_level": "cell",
                    "budget": budget,
                    "block_size_bp": block_size,
                    "exact_cell_precision": exact_point["precision"],
                    "all_wrong_cell_mean_precision": float(np.mean([row["precision"] for row in wrong_points])),
                    "exact_minus_all_wrong_mean": float(exact_point["precision"] - np.mean([row["precision"] for row in wrong_points])),
                    "ci95_low": low,
                    "ci95_high": high,
                    "exact_cell_hits": exact_point["emvar_hits"],
                    "all_wrong_cell_mean_hits": float(np.mean([row["emvar_hits"] for row in wrong_points])),
                    "wrong_prediction_contexts": 4,
                    "bootstrap_replicates": args.bootstrap_replicates,
                    "bootstrap_seed": bootstrap_seed,
                })
            macro_difference = np.mean(wrong_cell_differences, axis=0)
            low, high = interval(macro_difference)
            exact_points = [points[("context", cell, budget, cell)] for cell in CONFIRMATION_CELLS]
            wrong_point_means = []
            for outcome_cell in CONFIRMATION_CELLS:
                wrong_point_means.append(np.mean([
                    points[("context", outcome_cell, budget, prediction_cell)]["precision"]
                    for prediction_cell in CELLS
                    if prediction_cell != outcome_cell
                ]))
            wrong_contrasts.append({
                "post_hoc": "true",
                "cell_context": "confirmation_macro",
                "summary_level": "unweighted_four_cell_macro",
                "budget": budget,
                "block_size_bp": block_size,
                "exact_cell_precision": float(np.mean([row["precision"] for row in exact_points])),
                "all_wrong_cell_mean_precision": float(np.mean(wrong_point_means)),
                "exact_minus_all_wrong_mean": float(np.mean([row["precision"] for row in exact_points]) - np.mean(wrong_point_means)),
                "ci95_low": low,
                "ci95_high": high,
                "exact_cell_hits": sum(row["emvar_hits"] for row in exact_points),
                "all_wrong_cell_mean_hits": float(np.mean(wrong_point_means) * budget * len(CONFIRMATION_CELLS)),
                "wrong_prediction_contexts": 16,
                "bootstrap_replicates": args.bootstrap_replicates,
                "bootstrap_seed": bootstrap_seed,
            })

            macro_draws = {}
            for focal, reference in single_comparisons:
                macro_draws[f"{focal}_minus_{reference}"] = np.mean([
                    draws[("single", cell, budget, focal)] - draws[("single", cell, budget, reference)]
                    for cell in CONFIRMATION_CELLS
                ], axis=0)
            for replicate in range(args.bootstrap_replicates):
                bootstrap_records.append({
                    "post_hoc": "true",
                    "block_size_bp": block_size,
                    "bootstrap_seed": bootstrap_seed,
                    "replicate": replicate + 1,
                    "budget": budget,
                    **{name: values[replicate] for name, values in macro_draws.items()},
                    "exact_minus_all_wrong_mean": macro_difference[replicate],
                    **{
                        f"{cell}_exact_minus_all_wrong_mean": wrong_cell_differences[index][replicate]
                        for index, cell in enumerate(CONFIRMATION_CELLS)
                    },
                })

            block_label = "one_mb" if block_size == 1_000_000 else "five_mb"
            macro_estimate_lookup = {
                row["score"]: row
                for row in single_estimates
                if row["cell_context"] == "confirmation_macro"
                and row["budget"] == budget
                and row["block_size_bp"] == block_size
            }
            macro_contrast_lookup = {
                (row["score"], row["reference_score"]): row
                for row in single_contrasts
                if row["cell_context"] == "confirmation_macro"
                and row["budget"] == budget
                and row["block_size_bp"] == block_size
            }
            wrong_macro = next(
                row for row in wrong_contrasts
                if row["cell_context"] == "confirmation_macro"
                and row["budget"] == budget
                and row["block_size_bp"] == block_size
            )
            headline[str(budget)][block_label] = {
                "macro_precision_by_score": {
                    score: macro_estimate_lookup[score]["precision"] for score in SCORES
                },
                "macro_hits_by_score": {
                    score: macro_estimate_lookup[score]["emvar_hits"] for score in SCORES
                },
                "single_track_contrasts": {
                    f"{focal}_minus_{reference}": {
                        "precision_difference": macro_contrast_lookup[(focal, reference)]["precision_difference"],
                        "ci95": [macro_contrast_lookup[(focal, reference)]["ci95_low"], macro_contrast_lookup[(focal, reference)]["ci95_high"]],
                    }
                    for focal, reference in single_comparisons
                },
                "exact_cell_minus_all_wrong_mean": {
                    "exact_cell_precision": wrong_macro["exact_cell_precision"],
                    "all_wrong_cell_mean_precision": wrong_macro["all_wrong_cell_mean_precision"],
                    "precision_difference": wrong_macro["exact_minus_all_wrong_mean"],
                    "ci95": [wrong_macro["ci95_low"], wrong_macro["ci95_high"]],
                },
            }

    write_tsv(output / "single_track_estimates.tsv", single_estimates)
    write_tsv(output / "single_track_contrasts.tsv", single_contrasts)
    write_tsv(output / "wrong_cell_grid.tsv", wrong_grid)
    write_tsv(output / "wrong_cell_contrasts.tsv", wrong_contrasts)
    write_tsv_gz_deterministic(output / "bootstrap_macro_contrasts.tsv.gz", bootstrap_records)

    table_paths = sorted(path for path in output.iterdir() if path.is_file())
    summary = {
        "analysis_status": "post_hoc",
        "post_hoc_request_date": POST_HOC_DATE,
        "original_confirmatory_record_unchanged": True,
        "manuscript_unchanged": True,
        "original_commit": ORIGINAL_COMMIT,
        "confirmation_cells": list(CONFIRMATION_CELLS),
        "prediction_contexts": list(CELLS),
        "a549_role": "prediction_context_only; not a confirmatory outcome cell",
        "fixed_budgets": list(BUDGETS),
        "block_sizes_bp": list(BLOCK_SIZES),
        "bootstrap": {
            "replicates": args.bootstrap_replicates,
            "base_seed": SEED,
            "seed_by_block_size": {str(size): SEED + size for size in BLOCK_SIZES},
            "method": "chromosome-stratified block bootstrap conditional on each fixed ranked selection",
            "replicate_count_basis": "existing repository practice uses 2,000 replicates",
        },
        "single_track_reconstruction": reconstruction,
        "eligibility": eligibility,
        "authoritative_combined_and_avi_replay": {
            "checks": authoritative_checks,
            "all_passed": all(row["passed"] for row in authoritative_checks),
        },
        "headline": headline,
        "input_sha256": {
            name: sha256(path)
            for name, path in paths.items()
            if name not in ("utility",)
        } | {"authoritative_budget_utility": sha256(paths["utility"])},
        "analysis_script_sha256": sha256(Path(__file__)),
        "output_sha256": {path.name: sha256(path) for path in table_paths},
        "coverage_checks": {
            "single_track_estimate_rows": len(single_estimates),
            "single_track_contrast_rows": len(single_contrasts),
            "wrong_cell_grid_rows": len(wrong_grid),
            "wrong_cell_contrast_rows": len(wrong_contrasts),
            "bootstrap_draw_rows": len(bootstrap_records),
            "full_four_outcome_by_five_prediction_context_grid": True,
            "all_estimates_finite": all(
                math.isfinite(float(row[field]))
                for records, fields in (
                    (single_estimates, ("precision", "ci95_low", "ci95_high")),
                    (single_contrasts, ("precision_difference", "ci95_low", "ci95_high")),
                    (wrong_grid, ("precision", "ci95_low", "ci95_high")),
                    (wrong_contrasts, ("exact_minus_all_wrong_mean", "ci95_low", "ci95_high")),
                )
                for row in records
                for field in fields
            ),
        },
    }
    (output / "context_control_summary.json").write_text(
        json.dumps(json_ready(summary), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "post_hoc": True,
        "bootstrap_replicates": args.bootstrap_replicates,
        "headline": headline,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
