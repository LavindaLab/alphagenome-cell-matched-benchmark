#!/usr/bin/env python3
"""Independent structural and numerical verification of manuscript analysis outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path


CELLS = {"HepG2", "K562", "HCT116", "SK-N-SH"}
BUDGETS = {96, 384, 1000, 5000}
BLOCK_SIZES = {1_000_000, 5_000_000}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text())


def load_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(left, right, tolerance=1e-12):
    return math.isclose(float(left), float(right), rel_tol=0, abs_tol=tolerance)


def main() -> None:
    project = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    analysis = project / "manuscript/analysis"
    checks = []

    def check(name, condition, detail=""):
        checks.append({"check": name, "passed": bool(condition), "detail": detail})

    manifest = load_json(analysis / "manifest.json")
    for name, expected in manifest["output_sha256"].items():
        check(f"hash_{name}", sha256(analysis / name) == expected)
    check("analysis_script_bound", sha256(project / "tools/analyze_manuscript.py") == manifest["analysis_script_sha256"])
    check("bootstrap_replicates", manifest["bootstrap_replicates"] == 2000)
    check("scientific_rules_unchanged", manifest["scientific_rules_mutated"] is False)
    check("gate_outcomes_excluded", manifest["gate_outcomes_used"] is False)

    summary = load_json(analysis / "manuscript_analysis_summary.json")
    frozen = load_json(project / "confirmation/results/confirmation_results.json")
    check("frozen_primary_point_preserved", close(summary["frozen_primary"]["primary_macro_difference"], frozen["primary_macro_difference"]))
    check("frozen_primary_1mb_preserved", summary["frozen_primary"]["primary_1mb_bootstrap_95_ci"] == frozen["primary_1mb_bootstrap_95_ci"])
    check("frozen_primary_5mb_preserved", summary["frozen_primary"]["primary_5mb_bootstrap_95_ci"] == frozen["primary_5mb_bootstrap_95_ci"])
    check("failed_concordance_preserved", summary["frozen_primary"]["quantitative_value_screens"]["B_concordance_policy"] is False)
    check("primary_1mb_contains_zero", frozen["primary_1mb_bootstrap_95_ci"][0] <= 0 <= frozen["primary_1mb_bootstrap_95_ci"][1])
    check("primary_5mb_contains_zero", frozen["primary_5mb_bootstrap_95_ci"][0] <= 0 <= frozen["primary_5mb_bootstrap_95_ci"][1])
    check("no_new_threshold_optimization", summary["new_threshold_optimized_on_outcomes"] is False)

    frozen_utility = load_tsv(project / "confirmation/results/budget_utility.tsv")
    utility = {(row["cell_context"], int(row["budget"]), row["policy"]): row for row in frozen_utility}
    contrasts = load_tsv(analysis / "budget_precision_contrasts.tsv")
    central = [row for row in contrasts if row["comparator"] == "avi"]
    for row in central:
        if row["cell_context"] == "confirmation_macro":
            continue
        cell, budget = row["cell_context"], int(row["budget"])
        expected = float(utility[(cell, budget, "magnitude")]["precision"]) - float(utility[(cell, budget, "avi")]["precision"])
        check(f"central_point_{cell}_{budget}_{row['block_size']}", close(row["precision_difference"], expected))
    central_confirmatory = [
        row for row in central
        if row["cell_context"] in CELLS and int(row["budget"]) in BUDGETS and int(row["block_size"]) in BLOCK_SIZES
    ]
    check("central_complete_4x4x2", len(central_confirmatory) == 32)
    check("central_all_point_directions_positive", all(float(row["precision_difference"]) > 0 for row in central_confirmatory))
    check("central_all_intervals_positive", all(float(row["ci95_low"]) > 0 for row in central_confirmatory))
    macro_central = [row for row in central if row["cell_context"] == "confirmation_macro"]
    check("central_macro_complete_4x2", len(macro_central) == 8)
    check("central_macro_all_intervals_positive", all(float(row["ci95_low"]) > 0 for row in macro_central))

    for budget, expected_difference, expected_hits in ((384, 0.271484375, 417), (1000, 0.24325, 973)):
        row = next(
            item for item in macro_central
            if int(item["budget"]) == budget and int(item["block_size"]) == 1_000_000
        )
        check(f"headline_difference_{budget}", close(row["precision_difference"], expected_difference))
        check(f"headline_hits_{budget}", int(row["additional_hits"]) == expected_hits)

    baseline_macro = [row for row in contrasts if row["cell_context"] == "confirmation_macro" and int(row["block_size"]) == 1_000_000]
    check("frozen_baselines_complete", {(int(row["budget"]), row["comparator"]) for row in baseline_macro} == {
        (budget, comparator) for budget in BUDGETS for comparator in {"random", "pip", "avi", "concordance", "concordance_block_cap"}
    })

    bootstrap_draws = load_tsv(analysis / "budget_precision_bootstrap_draws.tsv.gz")
    check("budget_draw_count", len(bootstrap_draws) == 2 * 2000 * 4)
    check("budget_draw_dimensions", {(int(row["block_size"]), int(row["budget"])) for row in bootstrap_draws} == {
        (size, budget) for size in BLOCK_SIZES for budget in BUDGETS
    })

    rank_contrasts = load_tsv(analysis / "avi_context_rank_contrasts.tsv")
    rank_confirmatory = [row for row in rank_contrasts if row["cell_context"] in CELLS]
    check("rank_discordance_complete_4x2", len(rank_confirmatory) == 8)
    check("rank_discordance_all_cells_positive", all(float(row["prevalence_difference"]) > 0 for row in rank_confirmatory))
    check("rank_discordance_all_intervals_positive", all(float(row["difference_ci95_low"]) > 0 for row in rank_confirmatory))
    check("rank_discordance_all_ratios_above_one", all(float(row["ratio_ci95_low"]) > 1 for row in rank_confirmatory))

    indicators = load_tsv(analysis / "uncertainty_indicator_models.tsv")
    confirmatory_indicators = [row for row in indicators if row["confirmation_role"] == "confirmatory"]
    check("indicator_fdr_family_20", len(confirmatory_indicators) == 20)
    check("indicator_q_values_valid", all(0 <= float(row["q_value_bh_20_tests"]) <= 1 for row in confirmatory_indicators))
    rank_indicator = [row for row in confirmatory_indicators if row["indicator"] == "magnitude_only_vs_avi_only_top10"]
    check("rank_indicator_four_cells", len(rank_indicator) == 4)
    check("rank_indicator_reproducible", all(float(row["ci95_low"]) > 1 and float(row["q_value_bh_20_tests"]) < 0.05 for row in rank_indicator))

    ranking = load_tsv(analysis / "ranking_metrics.tsv")
    for cell in CELLS:
        values = {row["score"]: row for row in ranking if row["cell_context"] == cell}
        check(f"ranking_metrics_{cell}", set(values) == {"magnitude", "avi", "pip"})
        check(f"magnitude_ap_above_avi_{cell}", float(values["magnitude"]["average_precision"]) > float(values["avi"]["average_precision"]))
        check(f"magnitude_auc_above_avi_{cell}", float(values["magnitude"]["auroc"]) > float(values["avi"]["auroc"]))

    integrity = load_json(analysis / "integrity_audit.json")
    check("integrity_all_passed", integrity["all_checks_passed"] is True)
    check("protocol_hash_match", integrity["protocol_sha256"] == integrity["protocol_sha256_authorized"])
    check("evidence_hash_match", integrity["evidence_ledger_sha256"] == integrity["evidence_ledger_sha256_authorized"])
    check("no_duplicate_variant_cells", integrity["duplicate_variant_cell_rows"] == 0)
    check("atlas_complete", integrity["atlas_complete"] == 246674 and integrity["atlas_incomplete"] == 0)
    check("prior_result_checks", integrity["result_checks"] == "25/25")
    check("prior_atlas_checks", integrity["atlas_checks"] == "273/273")

    failed = [row for row in checks if not row["passed"]]
    receipt = {
        "passed": not failed,
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "failed_checks": failed,
        "frozen_primary_authoritative": True,
        "failed_concordance_hypothesis_preserved": True,
        "scientific_rules_mutated": False,
        "gate_outcomes_used": False,
    }
    path = analysis / "inspection_receipt.json"
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"passed={receipt['passed']} checks={receipt['checks_passed']}/{receipt['checks_total']}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
