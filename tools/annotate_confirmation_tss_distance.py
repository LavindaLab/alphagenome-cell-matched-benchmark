#!/usr/bin/env python3
"""Annotate frozen confirmation SNVs by nearest protein-coding GENCODE TSS."""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ATTR = re.compile(r'(\w+) "([^"]*)"')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest_gzip", type=Path)
    parser.add_argument("gencode_gtf_gzip", type=Path)
    parser.add_argument("output_gzip", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    tss = collections.defaultdict(set)
    with gzip.open(args.gencode_gtf_gzip, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = dict(ATTR.findall(fields[8]))
            if attrs.get("transcript_type") != "protein_coding":
                continue
            position = int(fields[3]) if fields[6] == "+" else int(fields[4])
            tss[fields[0]].add(position)
    sorted_tss = {chrom: sorted(values) for chrom, values in tss.items()}
    counts = collections.Counter()
    records = 0
    with gzip.open(args.manifest_gzip, "rt", newline="") as source, gzip.open(args.output_gzip, "wt", newline="") as output:
        reader = csv.DictReader(source, delimiter="\t")
        fields = ["variant_id_hg38", "nearest_protein_coding_tss_distance_bp", "promoter_proximal_2kb"]
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            positions = sorted_tss[row["chromosome"]]
            variant_position = int(row["position_1based"])
            index = bisect.bisect_left(positions, variant_position)
            candidates = positions[max(0, index - 1):min(len(positions), index + 1)]
            distance = min(abs(variant_position - value) for value in candidates)
            proximal = distance <= 2000
            writer.writerow({
                "variant_id_hg38": row["variant_id"],
                "nearest_protein_coding_tss_distance_bp": distance,
                "promoter_proximal_2kb": str(proximal).lower(),
            })
            counts["promoter_proximal" if proximal else "distal"] += 1
            records += 1
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "gencode_release": "50",
        "gencode_sha256": sha256(args.gencode_gtf_gzip),
        "manifest_sha256": sha256(args.manifest_gzip),
        "variants": records,
        "classification_counts": dict(sorted(counts.items())),
        "definition": "promoter-proximal if within 2,000 bp of any GENCODE 50 protein-coding transcript TSS; distal otherwise",
        "purpose": "prespecified secondary experimental-budget balance only",
        "output_sha256": sha256(args.output_gzip),
        "confirmation_mpra_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"variants={records} promoter={counts['promoter_proximal']} distal={counts['distal']}")


if __name__ == "__main__":
    main()
