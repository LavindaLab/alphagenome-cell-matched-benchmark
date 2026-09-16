#!/usr/bin/env python3
"""Create an outcome-independent retry manifest from incomplete API statuses."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("full_manifest", type=Path)
    parser.add_argument("status_tsv", type=Path)
    parser.add_argument("retry_manifest", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    with args.full_manifest.open(newline="") as handle:
        full = list(csv.DictReader(handle, delimiter="\t"))
    with args.status_tsv.open(newline="") as handle:
        status = list(csv.DictReader(handle, delimiter="\t"))
    incomplete = {row["variant_id_hg38"] for row in status if row["status"] == "incomplete"}
    retry = [row for row in full if row["variant_id"] in incomplete]
    if len(retry) != len(incomplete):
        raise ValueError("incomplete status IDs do not map one-to-one to manifest")
    with args.retry_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(full[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(retry)
    payload = args.retry_manifest.read_bytes()
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": "all and only incomplete IDs from prior outcome-blind API status",
        "retry_variants": len(retry),
        "retry_manifest_sha256": hashlib.sha256(payload).hexdigest(),
        "mpra_outcomes_used_for_selection": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"retry_variants={len(retry)}")


if __name__ == "__main__":
    main()
