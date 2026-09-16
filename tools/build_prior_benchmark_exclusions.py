#!/usr/bin/env python3
"""Construct exact confirmation exclusions from frozen prior-work coordinates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_hg38", type=Path)
    parser.add_argument("intervals_tsv", type=Path)
    parser.add_argument("explicit_variants_tsv", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    intervals = read(args.intervals_tsv)
    explicit = {row["variant_id_hg38"]: row for row in read(args.explicit_variants_tsv)}
    output = {}
    for row in read(args.candidate_hg38):
        variant = row["variant_id"]
        if variant in explicit:
            output[variant] = {
                "variant_id_hg38": variant,
                "source": explicit[variant]["source"],
                "detail": explicit[variant]["detail"],
            }
        position = int(row["position_1based"])
        for interval in intervals:
            if (
                row["chromosome"] == interval["chromosome"]
                and int(interval["start_1based"]) <= position <= int(interval["end_1based"])
            ):
                output[variant] = {
                    "variant_id_hg38": variant,
                    "source": interval["source"],
                    "detail": f"{interval['case']}:{interval['chromosome']}:{interval['start_1based']}-{interval['end_1based']}",
                }
    records = sorted(output.values(), key=lambda row: row["variant_id_hg38"])
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["variant_id_hg38", "source", "detail"]
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "exact pre-query exclusion of prior AlphaGenome/Atlas examples and CAGI5 intervals",
        "candidate_variants": sum(1 for _ in read(args.candidate_hg38)),
        "frozen_intervals": len(intervals),
        "frozen_explicit_variants": len(explicit),
        "exact_candidate_exclusions": len(records),
        "interval_source": "official Kircher saturation-MPRA portal GRCh38 coordinates; AlphaGenome Extended Data Fig. 8 locus list",
        "intervals_sha256": hashlib.sha256(args.intervals_tsv.read_bytes()).hexdigest(),
        "explicit_variants_sha256": hashlib.sha256(args.explicit_variants_tsv.read_bytes()).hexdigest(),
        "output_sha256": hashlib.sha256(args.output_tsv.read_bytes()).hexdigest(),
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"exact_prior_benchmark_exclusions={len(records)}")


if __name__ == "__main__":
    main()
