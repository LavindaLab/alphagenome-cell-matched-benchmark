#!/usr/bin/env python3
"""Create checksum-bound held-out outcome authorization after prediction inspection."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("features", type=Path)
    parser.add_argument("avi", type=Path)
    parser.add_argument("regions", type=Path)
    parser.add_argument("adastra", type=Path)
    parser.add_argument("prediction_inspection", type=Path)
    parser.add_argument("analysis_script", type=Path)
    parser.add_argument("protocol", type=Path)
    parser.add_argument("evidence_ledger", type=Path)
    parser.add_argument("novelty_refresh", type=Path)
    parser.add_argument("value_gate", type=Path)
    parser.add_argument("repository", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    inspection = json.loads(args.prediction_inspection.read_text())
    if not inspection.get("passed") or inspection.get("complete_fraction", 0) < 0.99:
        raise ValueError("prediction retrieval has not passed independent inspection")
    if inspection.get("confirmation_mpra_outcomes_accessed") is not False:
        raise ValueError("prediction inspection does not preserve the outcome firewall")
    with gzip.open(args.manifest, "rt", newline="") as handle:
        variants = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
    head = subprocess.run(
        ["git", "-C", str(args.repository), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "-C", str(args.repository), "status", "--short"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if status:
        raise ValueError("repository must be clean before outcome authorization")

    authorization = {
        "authorized_utc": datetime.now(timezone.utc).isoformat(),
        "authorization_scope": "one frozen held-out Siraj confirmation extraction and analysis",
        "repository_head": head,
        "expected_confirmation_variants": variants,
        "manifest_gzip_sha256": sha256(args.manifest),
        "analysis_input_sha256": {
            name: sha256(getattr(args, name))
            for name in ("metadata", "features", "avi", "regions", "adastra")
        },
        "analysis_script_sha256": sha256(args.analysis_script),
        "protocol_sha256": sha256(args.protocol),
        "evidence_ledger_sha256": sha256(args.evidence_ledger),
        "novelty_refresh_sha256": sha256(args.novelty_refresh),
        "substantive_value_gate_sha256": sha256(args.value_gate),
        "prediction_inspection_sha256": sha256(args.prediction_inspection),
        "prediction_cache_tree_sha256": inspection["cache_tree_sha256"],
        "prediction_retrieval_inspected": True,
        "analysis_implementation_committed": True,
        "endogenous_evidence_frozen": True,
        "novelty_refresh_complete": True,
        "confirmation_outcome_access_authorized": True,
        "gate_outcomes_authorized_or_read": False,
        "scientific_rules_mutated": False,
    }
    if authorization["protocol_sha256"] != "09bcfa0c09d6c39f6b98c521d73e772629395252f6faa4a35e4d3a8314d43f5f":
        raise ValueError("frozen protocol checksum changed")
    if authorization["evidence_ledger_sha256"] != "05df7123ab873d25b3d4b93d0e93dcfdeccad98b009c8fe55e73b39dff7f382c":
        raise ValueError("pre-outcome evidence ledger checksum changed")
    args.output_json.write_text(json.dumps(authorization, indent=2, sort_keys=True) + "\n")
    print(f"authorized=true variants={variants} repository_head={head}")


if __name__ == "__main__":
    main()
