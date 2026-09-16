#!/usr/bin/env python3
"""Managed completion of the frozen held-out evaluation after outcome-blind inspection."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_TABLE4_SHA256 = "723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2"
EXPECTED_BRANCH = "codex/alphagenome-academic-study"
PREDICTION_PATHS = (
    "research/alphagenome-study1/confirmation/atlas/retrieval_summary.json",
    "research/alphagenome-study1/confirmation/atlas/inspection_receipt.json",
    "research/alphagenome-study1/confirmation/prediction_features.tsv.gz",
    "research/alphagenome-study1/confirmation/avi_scores.tsv.gz",
    "research/alphagenome-study1/confirmation/atlas_variant_status.tsv.gz",
    "research/alphagenome-study1/confirmation/prediction_derivation_summary.json",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write(path: Path, value: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(value)
    temporary.replace(path)


def sha256(path: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repository", type=Path)
    parser.add_argument("postprocess_run", type=Path)
    parser.add_argument("managed_run", type=Path)
    parser.add_argument("table4", type=Path)
    args = parser.parse_args()
    study = args.repository / "research/alphagenome-study1"
    confirmation = study / "confirmation"
    atlas = confirmation / "atlas"
    args.managed_run.mkdir(parents=True, exist_ok=True)

    def stage(value: str) -> None:
        write(args.managed_run / "CURRENT_STAGE", value + "\n")
        write(args.managed_run / "HEARTBEAT", now() + "\n")

    def event(message: str, **details) -> None:
        with (args.managed_run / "events.jsonl").open("a") as handle:
            handle.write(json.dumps({"utc": now(), "message": message, **details}, sort_keys=True) + "\n")

    def command(name: str, values: list[str], cwd: Path | None = None) -> None:
        stage(name)
        with (args.managed_run / f"{name}.log").open("w") as log:
            subprocess.run(values, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)

    write(args.managed_run / "RUN_PID", f"{os.getpid()}\n")
    write(args.managed_run / "RUN_STATUS", "running\n")
    stage("waiting_for_independent_prediction_inspection")
    event("completion_pipeline_started", confirmation_outcomes_accessed=False)
    try:
        while True:
            status = (args.postprocess_run / "RUN_STATUS").read_text().strip()
            if status == "completed":
                break
            if status == "failed":
                raise RuntimeError("prediction postprocess failed")
            upstream = (args.postprocess_run / "CURRENT_STAGE").read_text().strip()
            stage(f"waiting_for_independent_prediction_inspection {upstream}")
            time.sleep(30)

        inspection = json.loads((atlas / "inspection_receipt.json").read_text())
        if not inspection.get("passed") or inspection.get("complete_fraction", 0) < 0.99:
            raise RuntimeError("prediction inspection did not pass")
        if inspection.get("confirmation_mpra_outcomes_accessed") is not False:
            raise RuntimeError("prediction inspection violated the outcome firewall")
        branch = subprocess.run(
            ["git", "branch", "--show-current"], cwd=args.repository,
            check=True, capture_output=True, text=True,
        ).stdout.strip()
        if branch != EXPECTED_BRANCH:
            raise RuntimeError(f"unexpected branch: {branch}")
        status_lines = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"], cwd=args.repository,
            check=True, capture_output=True, text=True,
        ).stdout.splitlines()
        expected_status = {f"?? {path}" for path in PREDICTION_PATHS}
        if set(status_lines) != expected_status:
            raise RuntimeError(f"unexpected pre-authorization worktree state: {status_lines}")

        command("commit_prediction_freeze", ["git", "add", *PREDICTION_PATHS], args.repository)
        command("commit_prediction_freeze_commit", ["git", "commit", "-m", "research: freeze inspected AlphaGenome confirmation predictions"], args.repository)
        if subprocess.run(["git", "status", "--porcelain"], cwd=args.repository, check=True, capture_output=True, text=True).stdout:
            raise RuntimeError("repository not clean after prediction freeze")

        authorization = confirmation / "OUTCOME_ACCESS_AUTHORIZATION.json"
        command("authorize_confirmation_outcomes", [
            sys.executable,
            str(study / "tools/authorize_confirmation_outcomes.py"),
            str(confirmation / "variants_hg38.tsv.gz"),
            str(confirmation / "variant_metadata.tsv.gz"),
            str(confirmation / "prediction_features.tsv.gz"),
            str(confirmation / "avi_scores.tsv.gz"),
            str(confirmation / "tss_distance.tsv.gz"),
            str(study / "evidence/adastra_exact_cell_confirmation.tsv.gz"),
            str(atlas / "inspection_receipt.json"),
            str(study / "tools/analyze_confirmation.py"),
            str(study / "PROTOCOL_v0.1.md"),
            str(study / "EVIDENCE_LEDGER.md"),
            str(confirmation / "NOVELTY_ACCESS_GATE_REFRESH_2026-09-13.md"),
            str(confirmation / "SUBSTANTIVE_VALUE_GATE_v1.md"),
            str(args.repository),
            str(authorization),
        ])
        command("commit_outcome_authorization", ["git", "add", str(authorization.relative_to(args.repository))], args.repository)
        command("commit_outcome_authorization_commit", ["git", "commit", "-m", "research: authorize frozen AlphaGenome confirmation outcomes"], args.repository)

        if sha256(args.table4) != EXPECTED_TABLE4_SHA256:
            raise RuntimeError("Siraj Table 4 checksum changed")
        outcomes = confirmation / "mpra_outcomes_confirmation.tsv.gz"
        outcome_receipt = confirmation / "outcome_extraction_receipt.json"
        command("extract_held_out_outcomes", [
            sys.executable,
            str(study / "tools/extract_confirmation_outcomes.py"),
            str(args.table4),
            str(confirmation / "variants_hg38.tsv.gz"),
            str(authorization),
            str(outcomes),
            str(outcome_receipt),
        ])

        results = confirmation / "results"
        command("run_frozen_confirmation_analysis", [
            sys.executable,
            str(study / "tools/analyze_confirmation.py"),
            str(outcomes),
            str(confirmation / "variant_metadata.tsv.gz"),
            str(confirmation / "prediction_features.tsv.gz"),
            str(confirmation / "avi_scores.tsv.gz"),
            str(confirmation / "tss_distance.tsv.gz"),
            str(study / "evidence/adastra_exact_cell_confirmation.tsv.gz"),
            str(authorization),
            str(outcome_receipt),
            str(results),
        ])
        result_inspection = results / "inspection_receipt.json"
        command("independently_verify_confirmation_results", [
            sys.executable,
            str(study / "tools/verify_confirmation_results.py"),
            str(results),
            str(authorization),
            str(result_inspection),
        ])
        verified = json.loads(result_inspection.read_text())
        if not verified.get("passed"):
            raise RuntimeError("held-out result inspection failed")

        result_paths = [
            str(outcomes.relative_to(args.repository)),
            str(outcome_receipt.relative_to(args.repository)),
            str(results.relative_to(args.repository)),
        ]
        command("commit_verified_confirmation_results", ["git", "add", *result_paths], args.repository)
        command("commit_verified_confirmation_results_commit", ["git", "commit", "-m", "research: record frozen AlphaGenome confirmation results"], args.repository)
        if subprocess.run(["git", "status", "--porcelain"], cwd=args.repository, check=True, capture_output=True, text=True).stdout:
            raise RuntimeError("repository not clean after result commit")

        result = json.loads((results / "confirmation_results.json").read_text())
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=args.repository, check=True, capture_output=True, text=True).stdout.strip()
        run_inspection = {
            "inspected_utc": now(),
            "passed": True,
            "result_inspection_sha256": sha256(result_inspection),
            "result_checks": f"{verified['checks_passed']}/{verified['checks_total']}",
            "result_commit": head,
            "scientific_rules_mutated": False,
            "manuscript_expansion_performed": False,
        }
        write(args.managed_run / "INSPECTION_RECEIPT.json", json.dumps(run_inspection, indent=2, sort_keys=True) + "\n")
        report = (
            "# Frozen AlphaGenome confirmation evaluation\n\n"
            f"- Result commit: `{head}`\n"
            f"- Primary difference: {result['primary_macro_difference']:.6f}\n"
            f"- 1 Mb 95% CI: {result['primary_1mb_bootstrap_95_ci']}\n"
            f"- Quantitative value screens: `{result['quantitative_value_screens']}`\n"
            f"- Independent result checks: {verified['checks_passed']}/{verified['checks_total']}\n"
            "- Scientific rules mutated: no\n"
            "- Manuscript expansion: not performed\n"
            "- Next and only gate: substantive scientific-value judgment.\n"
        )
        write(args.managed_run / "FINAL_REPORT.md", report)
        stage("held_out_evaluation_complete_pending_substantive_value_judgment")
        write(args.managed_run / "RUN_STATUS", "completed\n")
        write(args.managed_run / "COMPLETED", now() + "\n")
        event("completion_pipeline_completed", confirmation_outcomes_accessed=True, manuscript_expansion=False)
    except Exception as exc:
        message = " ".join(f"{type(exc).__name__}: {exc}".split())
        write(args.managed_run / "FINAL_REPORT.md", f"# Frozen AlphaGenome confirmation evaluation\n\nFAILED: {message}\nManuscript expansion was not performed.\n")
        stage("failed")
        write(args.managed_run / "RUN_STATUS", "failed\n")
        write(args.managed_run / "FAILED", now() + "\n")
        event("completion_pipeline_failed", error=message)
        raise


if __name__ == "__main__":
    main()
