#!/usr/bin/env python3
"""Independently verify frozen pilot result tables and key policy estimates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truth(value: str) -> bool:
    return value == "true"


def direction(value: str) -> int:
    number = float(value)
    return (number > 0) - (number < 0)


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance)


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def evaluate(selected_rows, rank, predicted_direction, concordant_only=False):
    requested = math.ceil(len(selected_rows) * 0.10)
    eligible = selected_rows
    if concordant_only:
        eligible = [row for row in selected_rows if row["concordance_eligible"]]
    picked = sorted(
        eligible,
        key=lambda row: (-float(row[rank]), row["tie_break_sha256"]),
    )[:requested]
    error = statistics.fmean(
        int(int(row[predicted_direction]) != row["outcome_direction"])
        for row in picked
    )
    return requested, len(picked), error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pilot_dir", type=Path)
    parser.add_argument("output_receipt", type=Path)
    args = parser.parse_args()
    pilot = args.pilot_dir
    result_dir = pilot / "results"

    outcomes = rows(pilot / "mpra_outcomes_gate.tsv")
    noncoding = {
        row["variant_id_hg38"]: truth(row["noncoding_eligible"])
        for row in rows(pilot / "noncoding_eligibility.tsv")
    }
    features = {
        (row["variant_id_hg38"], row["cell_context"]): row
        for row in rows(pilot / "prediction_features.tsv")
    }
    reconstructed = {cell: [] for cell in CELLS}
    for outcome in outcomes:
        key = (outcome["variant_id_hg38"], outcome["cell_context"])
        feature = features.get(key)
        if (
            feature is None
            or not noncoding[outcome["variant_id_hg38"]]
            or not truth(outcome["measured"])
            or not truth(outcome["enhancer_active"])
            or not truth(outcome["emvar"])
        ):
            continue
        reconstructed[outcome["cell_context"]].append(
            {
                "outcome_direction": direction(outcome["log2skew_alt_minus_ref"]),
                "magnitude_M": feature["magnitude_M"],
                "magnitude_predicted_direction": feature["magnitude_predicted_direction"],
                "concordance_C": feature["concordance_C"],
                "concordance_predicted_direction": feature["concordance_predicted_direction"],
                "concordance_eligible": truth(feature["concordance_eligible"]),
                "tie_break_sha256": feature["tie_break_sha256"],
            }
        )

    reported_rows = {
        row["cell_context"]: row for row in rows(result_dir / "primary_by_cell.tsv")
    }
    checks = []
    independent_differences = []
    for cell in CELLS:
        source = reconstructed[cell]
        requested_m, selected_m, error_m = evaluate(
            source, "magnitude_M", "magnitude_predicted_direction"
        )
        requested_c, selected_c, error_c = evaluate(
            source,
            "concordance_C",
            "concordance_predicted_direction",
            concordant_only=True,
        )
        reported = reported_rows[cell]
        independent_differences.append(error_c - error_m)
        checks.extend(
            [
                {"name": f"{cell}_population", "pass": len(source) == int(reported["population"])},
                {"name": f"{cell}_requested", "pass": requested_m == requested_c == int(reported["requested"])},
                {"name": f"{cell}_magnitude_selected", "pass": selected_m == int(reported["magnitude_selected"])},
                {"name": f"{cell}_concordance_selected", "pass": selected_c == int(reported["concordance_selected"])},
                {"name": f"{cell}_magnitude_error", "pass": close(error_m, float(reported["magnitude_error"]))},
                {"name": f"{cell}_concordance_error", "pass": close(error_c, float(reported["concordance_error"]))},
            ]
        )

    result = json.loads((result_dir / "pilot_results.json").read_text())
    independent_macro = statistics.fmean(independent_differences)
    checks.append({
        "name": "primary_macro_difference",
        "pass": close(independent_macro, float(result["primary_macro_difference"])),
    })

    expected_rows = {
        "primary_by_cell.tsv": 5,
        "risk_coverage.tsv": 25,
        "precision_lift.tsv": 50,
        "budget_utility.tsv": 120,
        "discordance.tsv": 10,
        "bootstrap_draws.tsv": 4000,
    }
    for filename, expected in expected_rows.items():
        checks.append({
            "name": f"{filename}_row_count",
            "pass": len(rows(result_dir / filename)) == expected,
        })

    draws = rows(result_dir / "bootstrap_draws.tsv")
    for block_size, result_key, draw_key in (
        (1_000_000, "primary_1mb_bootstrap_95_ci", "macro_difference"),
        (5_000_000, "primary_5mb_bootstrap_95_ci", "macro_difference"),
        (1_000_000, "continuation_1mb_bootstrap_95_ci", "continuation_macro_difference"),
        (5_000_000, "continuation_5mb_bootstrap_95_ci", "continuation_macro_difference"),
    ):
        values = [
            float(row[draw_key]) for row in draws
            if int(row["block_size"]) == block_size
        ]
        calculated = [percentile(values, 0.025), percentile(values, 0.975)]
        reported = result[result_key]
        checks.append({
            "name": f"{result_key}_from_draws",
            "pass": close(calculated[0], reported[0]) and close(calculated[1], reported[1]),
        })

    utility = rows(result_dir / "budget_utility.tsv")
    utility_checks = []
    for row in utility:
        selected = int(row["selected"])
        hits = int(row["emvar_hits"])
        utility_checks.append(
            selected <= min(int(row["budget"]), int(row["eligible_population"]))
            and hits <= selected
            and int(row["independent_blocks_with_hit"]) <= hits
        )
    checks.append({"name": "budget_utility_internal_bounds", "pass": all(utility_checks)})

    file_hashes = {
        str(path.relative_to(pilot)): sha256(path)
        for path in sorted(
            [
                pilot / "mpra_outcomes_gate.tsv",
                pilot / "variant_metadata_gate.tsv",
                pilot / "noncoding_eligibility.tsv",
                pilot / "prediction_features.tsv",
                pilot / "avi_merged" / "avi_scores.jsonl",
                *[result_dir / name for name in expected_rows],
                result_dir / "pilot_results.json",
            ]
        )
    }
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "independent frozen-pilot result verification",
        "confirmation_split_accessed": False,
        "independent_primary_macro_difference": independent_macro,
        "checks": checks,
        "checks_passed": sum(item["pass"] for item in checks),
        "checks_total": len(checks),
        "pass": all(item["pass"] for item in checks),
        "file_sha256": file_hashes,
    }
    args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        f"pass={receipt['pass']} checks={receipt['checks_passed']}/{receipt['checks_total']} "
        f"primary_macro_difference={independent_macro:.9f}"
    )
    raise SystemExit(0 if receipt["pass"] else 1)


if __name__ == "__main__":
    main()
