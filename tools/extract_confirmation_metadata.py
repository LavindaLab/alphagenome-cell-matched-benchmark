#!/usr/bin/env python3
"""Extract outcome-blind Table 1 metadata for frozen confirmation variants."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from extract_gate_outcomes import extract_table1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table1", type=Path)
    parser.add_argument("manifest_gzip", type=Path)
    parser.add_argument("output_gzip", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    with gzip.open(args.manifest_gzip, "rt", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if not manifest or any(row["frozen_split"] != "confirmation" for row in manifest):
        raise ValueError("manifest is not the frozen confirmation split")
    final_by_source = {row["source_variant_id_hg19"]: row for row in manifest}
    if len(final_by_source) != len(manifest):
        raise ValueError("duplicate source identifiers")
    metadata = extract_table1(args.table1, final_by_source)
    args.output_gzip.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output_gzip, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metadata)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind PIP/category baseline metadata",
        "table1_sha256": sha256(args.table1),
        "manifest_sha256": sha256(args.manifest_gzip),
        "variants": len(metadata),
        "fine_mapped_candidates": sum(row["fine_mapped_candidate"] == "true" for row in metadata),
        "pip_available": sum(row["pip_available"] == "true" for row in metadata),
        "fields_accessed": ["Variant", "PIP", "Category", "Variant Selection Type", "CS ID"],
        "table4_opened": False,
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
        "output_sha256": sha256(args.output_gzip),
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        f"metadata={len(metadata)} fine_mapped={receipt['fine_mapped_candidates']} "
        f"pip={receipt['pip_available']}"
    )


if __name__ == "__main__":
    main()
