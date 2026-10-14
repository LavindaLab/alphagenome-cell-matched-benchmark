#!/usr/bin/env python3
"""Independently verify manuscript claims, figure provenance and package boundaries."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"
CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(left, right, tolerance=1e-12):
    return math.isclose(float(left), float(right), rel_tol=0, abs_tol=tolerance)


def main() -> None:
    checks = []

    def check(name, condition):
        checks.append({"name": name, "passed": bool(condition)})

    required = [
        "MANUSCRIPT.md", "SUPPLEMENTARY_METHODS_RESULTS.md",
        "CLAIMS_EVIDENCE_TABLE.tsv", "LITERATURE_CLAIM_LOCK_2026-09-14.md",
        "figure_manifest.json",
    ]
    for name in required:
        check(f"required_{name}", (MANUSCRIPT / name).is_file())

    analysis_receipt = json.loads((MANUSCRIPT / "analysis/inspection_receipt.json").read_text())
    integrity = json.loads((MANUSCRIPT / "analysis/integrity_audit.json").read_text())
    primary = json.loads((PROJECT / "confirmation/results/confirmation_results.json").read_text())
    check("analysis_independent_100_of_100", analysis_receipt["passed"] and analysis_receipt["checks_passed"] == 100)
    check("integrity_all_pass", integrity["all_checks_passed"])
    check("gate_outcomes_not_used", not integrity["gate_outcomes_read"] and not analysis_receipt["gate_outcomes_used"])
    check("scientific_rules_unchanged", not integrity["scientific_rules_mutated"] and not analysis_receipt["scientific_rules_mutated"])
    check("atlas_complete", integrity["atlas_complete"] == 246_674 and integrity["atlas_incomplete"] == 0)
    check("no_duplicate_variant_cells", integrity["duplicate_variant_cell_rows"] == 0)
    check("primary_negative_preserved", analysis_receipt["failed_concordance_hypothesis_preserved"])
    check("primary_ci_includes_zero", primary["primary_1mb_bootstrap_95_ci"][0] <= 0 <= primary["primary_1mb_bootstrap_95_ci"][1])

    contrasts = rows(MANUSCRIPT / "analysis/budget_precision_contrasts.tsv")
    central = {
        int(row["budget"]): row for row in contrasts
        if row["cell_context"] == "confirmation_macro" and row["comparator"] == "avi" and int(row["block_size"]) == 1_000_000
    }
    check("budget_384_point", close(central[384]["precision_difference"], 0.271484375))
    check("budget_384_interval", close(central[384]["ci95_low"], 0.22934001956744113) and close(central[384]["ci95_high"], 0.31222111594203505))
    check("budget_1000_point", close(central[1000]["precision_difference"], 0.24324999999999997))
    check("budget_1000_interval", close(central[1000]["ci95_low"], 0.21867802087924407) and close(central[1000]["ci95_high"], 0.2682184419789291))
    for budget in (96, 384, 1000, 5000):
        for cell in CONFIRMATION_CELLS:
            row = next(item for item in contrasts if item["cell_context"] == cell and item["comparator"] == "avi" and int(item["budget"]) == budget and int(item["block_size"]) == 1_000_000)
            check(f"positive_cell_{cell}_{budget}", float(row["ci95_low"]) > 0)

    rank = next(
        row for row in rows(MANUSCRIPT / "analysis/avi_context_rank_contrasts.tsv")
        if row["cell_context"] == "confirmation_macro" and int(row["block_size"]) == 1_000_000
    )
    check("rank_disagreement_point", close(rank["prevalence_difference"], 0.08380094639236486))
    check("rank_disagreement_interval_positive", float(rank["difference_ci95_low"]) > 0)
    check("rank_disagreement_ratio", close(rank["prevalence_ratio"], 2.4123558685724773))

    utility = rows(PROJECT / "confirmation/results/budget_utility.tsv")
    totals = {}
    for budget in (384, 1000):
        for policy in ("magnitude", "avi"):
            selected = [row for row in utility if row["cell_context"] in CONFIRMATION_CELLS and int(row["budget"]) == budget and row["policy"] == policy]
            totals[(budget, policy)] = (
                sum(int(row["emvar_hits"]) for row in selected),
                sum(int(row["independent_blocks_with_hit"]) for row in selected),
            )
    check("yield_384", totals[(384, "magnitude")] == (995, 896) and totals[(384, "avi")] == (578, 529))
    check("yield_1000", totals[(1000, "magnitude")] == (2250, 1764) and totals[(1000, "avi")] == (1277, 1090))

    claims = rows(MANUSCRIPT / "CLAIMS_EVIDENCE_TABLE.tsv")
    check("claims_22", len(claims) == 22 and {row["claim_id"] for row in claims} == {f"C{i:02d}" for i in range(1, 23)})

    context_receipt = json.loads((MANUSCRIPT / "analysis/context-controls/verification_receipt.json").read_text())
    check("context_controls_15_of_15", context_receipt["passed"] and context_receipt["checks_passed"] == 15)
    context_rows = rows(MANUSCRIPT / "analysis/context-controls/wrong_cell_contrasts.tsv")
    context_macro = {
        int(row["budget"]): row for row in context_rows
        if row["cell_context"] == "confirmation_macro" and int(row["block_size_bp"]) == 1_000_000
    }
    check("context_384_point", close(context_macro[384]["exact_minus_all_wrong_mean"], 0.17545572916666674))
    check("context_1000_point", close(context_macro[1000]["exact_minus_all_wrong_mean"], 0.1500625))

    figure_manifest = json.loads((MANUSCRIPT / "figure_manifest.json").read_text())
    check("figure_count", figure_manifest["figure_count"] == 6 and figure_manifest["pdf_pages"] == 6)
    check("renderer_hash", figure_manifest["renderer_sha256"] == sha256(PROJECT / "tools/render_manuscript_figures.py"))
    for relative, expected in figure_manifest["files"].items():
        path = MANUSCRIPT / relative
        check(f"figure_hash_{relative}", path.is_file() and sha256(path) == expected)
    pdf_path = MANUSCRIPT / "output/pdf/AlphaGenome_Atlas_manuscript_figures.pdf"
    pdf_info = subprocess.run(["pdfinfo", str(pdf_path)], check=True, capture_output=True, text=True).stdout
    check("pdf_six_pages", re.search(r"^Pages:\s+6$", pdf_info, re.MULTILINE) is not None)
    check("six_svg", len(list((MANUSCRIPT / "figures").glob("*.svg"))) == 6)
    check("six_png", len(list((MANUSCRIPT / "figures").glob("*.png"))) == 6)
    check("thirteen_figure_source_tables", len(list((MANUSCRIPT / "source-data").glob("*.tsv"))) == 13)

    manuscript_text = (MANUSCRIPT / "MANUSCRIPT.md").read_text()
    for exact in ("0.1755", "0.1501", "417 additional", "973", "hypothesis failed", "post hoc"):
        check(f"manuscript_contains_{exact}", exact in manuscript_text)
    inspectable_suffixes = {".md", ".tsv", ".json", ".py", ".svg"}
    inspectable = [path for path in MANUSCRIPT.rglob("*") if path.is_file() and path.suffix in inspectable_suffixes]
    inspectable.extend((PROJECT / "tools" / name for name in ("render_manuscript_figures.py", "verify_manuscript_package.py")))
    searchable_text = "\n".join(path.read_text(errors="ignore") for path in inspectable)
    check("no_api_key_pattern", re.search(r"AIza[0-9A-Za-z_-]{30,}", searchable_text) is None)
    raw_atlas_patterns = ("cache", "prediction_features", "avi_scores", "atlas_variant_status")
    raw_atlas_files = [
        str(path.relative_to(MANUSCRIPT)) for path in MANUSCRIPT.rglob("*")
        if path.is_file() and any(pattern in path.name.lower() for pattern in raw_atlas_patterns)
    ]
    check("no_raw_atlas_outputs", not raw_atlas_files)

    failed = [item["name"] for item in checks if not item["passed"]]
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "passed": not failed,
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "failed_checks": failed,
        "manual_gates_remaining": [
            "posting-day novelty and terms refresh",
            "human citation and final submission-PDF approval",
            "public repository archive DOI before journal submission",
        ],
        "scientific_scope_complete": not failed,
        "primary_negative_preserved": analysis_receipt["failed_concordance_hypothesis_preserved"],
        "raw_atlas_outputs_in_manuscript_tree": raw_atlas_files,
    }
    (MANUSCRIPT / "PACKAGE_INSPECTION_RECEIPT.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"passed={receipt['passed']} checks={receipt['checks_passed']}/{receipt['checks_total']} failed={len(failed)}")
    if failed:
        raise SystemExit("failed checks: " + ", ".join(failed))


if __name__ == "__main__":
    main()
