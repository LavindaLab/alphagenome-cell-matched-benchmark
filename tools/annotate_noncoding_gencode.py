#!/usr/bin/env python3
"""Apply the frozen GENCODE protein-coding CDS/essential-splice exclusion."""

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
    parser.add_argument("variants_hg38_tsv", type=Path)
    parser.add_argument("gencode_gtf_gz", type=Path)
    parser.add_argument("eligibility_tsv", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    with args.variants_hg38_tsv.open(newline="") as handle:
        variants = list(csv.DictReader(handle, delimiter="\t"))
    by_chrom = collections.defaultdict(list)
    lookup = {}
    for row in variants:
        position = int(row["position_1based"])
        by_chrom[row["chromosome"]].append(position)
        lookup[(row["chromosome"], position)] = row["variant_id"]
    for positions in by_chrom.values():
        positions.sort()
    cds = set()
    exon_max = collections.defaultdict(int)
    candidate_boundaries = []
    with gzip.open(args.gencode_gtf_gz, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"CDS", "exon"}:
                continue
            attrs = dict(ATTR.findall(fields[8]))
            if attrs.get("transcript_type") != "protein_coding":
                continue
            chromosome = fields[0]
            positions = by_chrom.get(chromosome)
            if not positions:
                continue
            start, end = int(fields[3]), int(fields[4])
            if fields[2] == "CDS":
                left = bisect.bisect_left(positions, start)
                right = bisect.bisect_right(positions, end)
                cds.update((chromosome, position) for position in positions[left:right])
                continue
            transcript = attrs.get("transcript_id")
            exon_number = int(attrs.get("exon_number", "0"))
            if not transcript or not exon_number:
                continue
            exon_max[transcript] = max(exon_max[transcript], exon_number)
            strand = fields[6]
            for boundary, candidate_positions in (
                ("start", (start - 2, start - 1)),
                ("end", (end + 1, end + 2)),
            ):
                for position in candidate_positions:
                    if (chromosome, position) in lookup:
                        candidate_boundaries.append((chromosome, position, transcript, exon_number, strand, boundary))
    splice = set()
    for chromosome, position, transcript, exon_number, strand, boundary in candidate_boundaries:
        maximum = exon_max[transcript]
        internal = (
            (strand == "+" and boundary == "start" and exon_number > 1)
            or (strand == "+" and boundary == "end" and exon_number < maximum)
            or (strand == "-" and boundary == "start" and exon_number < maximum)
            or (strand == "-" and boundary == "end" and exon_number > 1)
        )
        if internal:
            splice.add((chromosome, position))
    output = []
    for row in variants:
        key = (row["chromosome"], int(row["position_1based"]))
        coding = key in cds
        essential_splice = key in splice
        output.append(
            {
                "variant_id_hg38": row["variant_id"],
                "protein_coding_cds_overlap": str(coding).lower(),
                "essential_splice_overlap": str(essential_splice).lower(),
                "noncoding_eligible": str(not coding and not essential_splice).lower(),
            }
        )
    with args.eligibility_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "gencode_source": "GENCODE release 50 comprehensive gene annotation GTF",
        "gencode_sha256": sha256(args.gencode_gtf_gz),
        "variants": len(output),
        "protein_coding_cds_overlap": sum(row["protein_coding_cds_overlap"] == "true" for row in output),
        "essential_splice_overlap": sum(row["essential_splice_overlap"] == "true" for row in output),
        "noncoding_eligible": sum(row["noncoding_eligible"] == "true" for row in output),
        "rule": "exclude protein-coding transcript CDS overlap or either of two intronic bases adjacent to an internal protein-coding exon boundary",
        "confirmation_split_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"variants={len(output)} noncoding={receipt['noncoding_eligible']} cds={receipt['protein_coding_cds_overlap']} splice={receipt['essential_splice_overlap']}")


if __name__ == "__main__":
    main()
