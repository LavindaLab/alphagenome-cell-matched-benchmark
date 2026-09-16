#!/usr/bin/env python3
"""Apply only technical, annotation, and frozen exact-overlap exclusions."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PRIMARY_CHROMOSOMES = {f"chr{number}" for number in range(1, 23)} | {"chrX", "chrY"}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, records) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def set_hash(records, field):
    payload = "\n".join(sorted(row[field] for row in records)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_hg19", type=Path)
    parser.add_argument("lifted_hg38", type=Path)
    parser.add_argument("reference_validation", type=Path)
    parser.add_argument("noncoding_eligibility", type=Path)
    parser.add_argument("exact_overlap_exclusions", type=Path)
    parser.add_argument("final_hg19", type=Path)
    parser.add_argument("final_hg38", type=Path)
    parser.add_argument("exclusions_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    args = parser.parse_args()

    source = {row["variant_id"]: row for row in read_tsv(args.source_hg19)}
    lifted = read_tsv(args.lifted_hg38)
    validation = {
        row["variant_id"]: row["reference_status"]
        for row in read_tsv(args.reference_validation)
    }
    noncoding = {
        row["variant_id_hg38"]: row["noncoding_eligible"] == "true"
        for row in read_tsv(args.noncoding_eligibility)
    }
    overlap_rows = read_tsv(args.exact_overlap_exclusions)
    overlap_ids = {row["variant_id_hg38"] for row in overlap_rows}
    by_target = collections.defaultdict(list)
    for row in lifted:
        by_target[row["variant_id"]].append(row)

    accepted = []
    exclusions = []
    lifted_sources = {row["source_variant_id_hg19"] for row in lifted}
    for source_id, source_row in source.items():
        if source_id not in lifted_sources:
            exclusions.append({"source_variant_id_hg19": source_id, "variant_id_hg38": "", "reason": "unmapped_or_nonunique_liftover"})
    for row in lifted:
        target = row["variant_id"]
        reason = None
        if row["chromosome"] not in PRIMARY_CHROMOSOMES:
            reason = "unsupported_nonprimary_hg38_contig"
        elif validation.get(target) != "match":
            reason = "hg38_reference_mismatch"
        elif len(by_target[target]) != 1:
            reason = "duplicate_hg38_allele_representation"
        elif not noncoding.get(target, False):
            reason = "protein_coding_CDS_or_essential_splice"
        elif target in overlap_ids:
            reason = "frozen_AlphaGenome_or_Atlas_exact_overlap"
        if reason:
            exclusions.append({"source_variant_id_hg19": row["source_variant_id_hg19"], "variant_id_hg38": target, "reason": reason})
        else:
            accepted.append(row)

    accepted.sort(key=lambda row: int(row["selection_rank"]))
    final_hg19 = [source[row["source_variant_id_hg19"]] for row in accepted]
    exclusions.sort(key=lambda row: (row["source_variant_id_hg19"], row["reason"]))
    write_tsv(args.final_hg19, final_hg19)
    write_tsv(args.final_hg38, accepted)
    write_tsv(args.exclusions_tsv, exclusions)
    counts = collections.Counter(row["reason"] for row in exclusions)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "final outcome-blind confirmation query manifest",
        "source_confirmation_snvs": len(source),
        "uniquely_lifted_rows": len(lifted),
        "final_query_variants": len(accepted),
        "exclusion_counts": dict(sorted(counts.items())),
        "primary_chromosome_rule": sorted(PRIMARY_CHROMOSOMES),
        "exact_overlap_list_sha256": hashlib.sha256(args.exact_overlap_exclusions.read_bytes()).hexdigest(),
        "final_hg19_set_sha256": set_hash(final_hg19, "variant_id"),
        "final_hg38_set_sha256": set_hash(accepted, "variant_id"),
        "all_reference_valid": all(validation[row["variant_id"]] == "match" for row in accepted),
        "all_noncoding_eligible": all(noncoding[row["variant_id"]] for row in accepted),
        "mpra_functional_outcomes_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"source={len(source)} final={len(accepted)} exclusions={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
