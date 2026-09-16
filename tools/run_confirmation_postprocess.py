#!/usr/bin/env python3
"""Managed outcome-blind derivation and inspection after confirmation retrieval."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write(path: Path, value: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(value)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repository", type=Path)
    parser.add_argument("retrieval_run", type=Path)
    parser.add_argument("managed_run", type=Path)
    args = parser.parse_args()
    study = args.repository / "research/alphagenome-study1"
    atlas = study / "confirmation/atlas"
    args.managed_run.mkdir(parents=True, exist_ok=True)

    def stage(value: str) -> None:
        write(args.managed_run / "CURRENT_STAGE", value + "\n")
        write(args.managed_run / "HEARTBEAT", now() + "\n")

    def event(message: str) -> None:
        with (args.managed_run / "events.jsonl").open("a") as handle:
            handle.write(json.dumps({"utc": now(), "message": message, "mpra_outcomes_accessed": False}, sort_keys=True) + "\n")

    write(args.managed_run / "RUN_PID", f"{os.getpid()}\n")
    write(args.managed_run / "RUN_STATUS", "running\n")
    stage("waiting_for_outcome_blind_retrieval")
    event("postprocess_started")
    try:
        while True:
            retrieval_status = (args.retrieval_run / "RUN_STATUS").read_text().strip()
            if retrieval_status == "completed":
                break
            if retrieval_status == "failed":
                raise RuntimeError("outcome-blind retrieval failed")
            retrieval_stage = (args.retrieval_run / "CURRENT_STAGE").read_text().strip()
            stage(f"waiting_for_outcome_blind_retrieval {retrieval_stage}")
            time.sleep(30)

        stage("deriving_outcome_blind_prediction_features")
        with (args.managed_run / "derive.log").open("w") as log:
            subprocess.run([
                sys.executable,
                str(study / "tools/derive_confirmation_prediction_features.py"),
                str(atlas),
                str(study / "confirmation/prediction_features.tsv.gz"),
                str(study / "confirmation/avi_scores.tsv.gz"),
                str(study / "confirmation/atlas_variant_status.tsv.gz"),
                str(study / "confirmation/prediction_derivation_summary.json"),
            ], check=True, stdout=log, stderr=subprocess.STDOUT)

        stage("independently_inspecting_prediction_cache_and_derivation")
        with (args.managed_run / "verify.log").open("w") as log:
            subprocess.run([
                sys.executable,
                str(study / "tools/verify_confirmation_retrieval.py"),
                str(study / "confirmation/variants_hg38.tsv.gz"),
                str(atlas),
                str(study / "confirmation/prediction_features.tsv.gz"),
                str(study / "confirmation/avi_scores.tsv.gz"),
                str(study / "confirmation/atlas_variant_status.tsv.gz"),
                str(study / "confirmation/prediction_derivation_summary.json"),
                str(atlas / "inspection_receipt.json"),
            ], check=True, stdout=log, stderr=subprocess.STDOUT)

        inspection = json.loads((atlas / "inspection_receipt.json").read_text())
        report = (
            "# AlphaGenome confirmation postprocess\n\n"
            f"- Independent inspection passed: {inspection['passed']}\n"
            f"- Checks: {inspection['checks_passed']}/{inspection['checks_total']}\n"
            f"- Complete variants: {inspection['complete']}/{inspection['manifest_variants']}\n"
            "- MPRA confirmation outcomes accessed: no\n"
            "- Next gate: commit inspected predictions, then generate checksum-bound outcome authorization.\n"
        )
        write(args.managed_run / "FINAL_REPORT.md", report)
        stage("prediction_derivation_and_independent_inspection_complete")
        write(args.managed_run / "RUN_STATUS", "completed\n")
        write(args.managed_run / "COMPLETED", now() + "\n")
        event("postprocess_completed")
    except Exception as exc:
        message = " ".join(f"{type(exc).__name__}: {exc}".split())
        write(args.managed_run / "FINAL_REPORT.md", f"# AlphaGenome confirmation postprocess\n\nFAILED: {message}\nMPRA outcomes were not accessed.\n")
        stage("failed")
        write(args.managed_run / "RUN_STATUS", "failed\n")
        write(args.managed_run / "FAILED", now() + "\n")
        event("postprocess_failed")
        raise


if __name__ == "__main__":
    main()
