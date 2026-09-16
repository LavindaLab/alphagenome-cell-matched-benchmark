#!/usr/bin/env python3
"""Independently verify frozen held-out confirmation result arithmetic."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def same_number(left, right, tolerance=1e-12):
    return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("authorization", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()

    result = json.loads((args.results_dir / "confirmation_results.json").read_text())
    manifest = json.loads((args.results_dir / "result_manifest.json").read_text())
    authorization = json.loads(args.authorization.read_text())
    checks = []

    def check(name, passed, observed=None, expected=None):
        checks.append({"check": name, "passed": bool(passed), "observed": observed, "expected": expected})

    for name, expected in manifest["output_sha256"].items():
        path = args.results_dir / name
        check(f"output_hash_{name}", path.is_file() and sha256(path) == expected, sha256(path) if path.is_file() else None, expected)
    check("authorization_hash", manifest["input_sha256"]["authorization"] == sha256(args.authorization))
    check("gate_outcomes_excluded", manifest["gate_outcomes_used"] is False)
    check("confirmation_outcomes_used", manifest["confirmation_outcomes_used"] is True)
    check("scientific_rules_unchanged", manifest["scientific_rules_mutated"] is False and result["scientific_rules_mutated"] is False)
    check("analysis_script_bound", authorization["analysis_script_sha256"] == manifest["input_sha256"]["analysis_script"])

    primary = read_tsv(args.results_dir / "primary_by_cell.tsv")
    primary_by_cell = {row["cell_context"]: row for row in primary}
    check("primary_five_cells", set(primary_by_cell) == {*CONFIRMATION_CELLS, "A549"})
    check("a549_descriptive", primary_by_cell["A549"]["confirmation_role"] == "descriptive_count_rule_downgrade")
    macro = statistics.fmean(float(primary_by_cell[cell]["difference"]) for cell in CONFIRMATION_CELLS)
    five_cell = statistics.fmean(float(row["difference"]) for row in primary)
    check("primary_macro", same_number(macro, result["primary_macro_difference"]), macro, result["primary_macro_difference"])
    check("five_cell_macro", same_number(five_cell, result["five_cell_descriptive_macro_difference"]), five_cell, result["five_cell_descriptive_macro_difference"])
    check("primary_counts", {cell: int(row["population"]) for cell, row in primary_by_cell.items()} == result["primary_counts_by_cell"])

    bootstrap = read_tsv(args.results_dir / "bootstrap_draws.tsv.gz")
    one_mb = [float(row["confirmation_macro_difference"]) for row in bootstrap if int(row["block_size"]) == 1_000_000]
    five_mb = [float(row["confirmation_macro_difference"]) for row in bootstrap if int(row["block_size"]) == 5_000_000]
    check("bootstrap_counts", len(one_mb) == len(five_mb) == 2000, [len(one_mb), len(five_mb)], [2000, 2000])
    ci_1mb = [float(np.quantile(one_mb, 0.025, method="linear")), float(np.quantile(one_mb, 0.975, method="linear"))]
    ci_5mb = [float(np.quantile(five_mb, 0.025, method="linear")), float(np.quantile(five_mb, 0.975, method="linear"))]
    check("bootstrap_1mb_ci", all(same_number(a, b) for a, b in zip(ci_1mb, result["primary_1mb_bootstrap_95_ci"], strict=True)), ci_1mb, result["primary_1mb_bootstrap_95_ci"])
    check("bootstrap_5mb_ci", all(same_number(a, b) for a, b in zip(ci_5mb, result["primary_5mb_bootstrap_95_ci"], strict=True)), ci_5mb, result["primary_5mb_bootstrap_95_ci"])

    utility = read_tsv(args.results_dir / "budget_utility.tsv")
    utility_lookup = {(row["cell_context"], int(row["budget"]), row["policy"]): row for row in utility}
    screens = {}
    for budget in (384, 1000):
        magnitude = [float(utility_lookup[(cell, budget, "magnitude")]["precision"]) for cell in CONFIRMATION_CELLS]
        avi = [float(utility_lookup[(cell, budget, "avi")]["precision"]) for cell in CONFIRMATION_CELLS]
        macro_magnitude, macro_avi = statistics.fmean(magnitude), statistics.fmean(avi)
        screens[str(budget)] = {
            "magnitude_macro_precision": macro_magnitude,
            "avi_macro_precision": macro_avi,
            "absolute_difference": macro_magnitude - macro_avi,
            "relative_gain": macro_magnitude / macro_avi - 1 if macro_avi else math.inf,
            "cells_favoring_magnitude": sum(left > right for left, right in zip(magnitude, avi, strict=True)),
        }
    for budget, observed in screens.items():
        expected = result["utility_screens"][budget]
        check(f"utility_screen_{budget}", all(
            observed[key] == expected[key] if isinstance(value, int) else same_number(observed[key], expected[key])
            for key, value in observed.items()
        ), observed, expected)
    screen_a = all(value["absolute_difference"] >= 0.03 and value["relative_gain"] >= 0.20 and value["cells_favoring_magnitude"] >= 3 for value in screens.values())

    primary_differences = [float(primary_by_cell[cell]["difference"]) for cell in CONFIRMATION_CELLS]
    same_direction_cells = sum((value > 0) == (macro > 0) and value != 0 for value in primary_differences)
    screen_b = abs(macro) >= 0.02 and (ci_1mb[0] > 0 or ci_1mb[1] < 0) and same_direction_cells >= 3

    discordance = read_tsv(args.results_dir / "discordance.tsv")
    discordance_lookup = {(row["cell_context"], row["dnase_h3k27ac_sign_disagreement"]): row for row in discordance}
    ratios = {}
    for cell in CONFIRMATION_CELLS:
        numerator = float(discordance_lookup[(cell, "false")]["emvar_prevalence"])
        denominator = float(discordance_lookup[(cell, "true")]["emvar_prevalence"])
        ratios[cell] = numerator / denominator if denominator > 0 else None
    check("agreement_ratios", all(
        observed is None and result["agreement_emvar_prevalence_ratio_by_confirmation_cell"][cell] is None
        or observed is not None and same_number(observed, result["agreement_emvar_prevalence_ratio_by_confirmation_cell"][cell])
        for cell, observed in ratios.items()
    ), ratios, result["agreement_emvar_prevalence_ratio_by_confirmation_cell"])
    screen_c = sum(value is not None and value >= 1.5 for value in ratios.values()) >= 3 and ci_1mb[0] > -0.02
    recomputed_screens = {
        "A_experimental_utility": screen_a,
        "B_concordance_policy": screen_b,
        "C_marker_policy_distinction": screen_c,
    }
    check("quantitative_value_screens", recomputed_screens == result["quantitative_value_screens"], recomputed_screens, result["quantitative_value_screens"])
    check("qualitative_gate_pending", result["qualitative_value_gate_pending"] is True)

    passed = all(item["passed"] for item in checks)
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "passed": passed,
        "checks_passed": sum(item["passed"] for item in checks),
        "checks_total": len(checks),
        "checks": checks,
        "primary_macro_difference": macro,
        "primary_1mb_bootstrap_95_ci": ci_1mb,
        "primary_5mb_bootstrap_95_ci": ci_5mb,
        "quantitative_value_screens": recomputed_screens,
        "gate_outcomes_used": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"passed={passed} checks={receipt['checks_passed']}/{receipt['checks_total']} screens={recomputed_screens}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
