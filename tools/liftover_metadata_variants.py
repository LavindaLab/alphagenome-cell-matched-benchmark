#!/usr/bin/env python3
"""Lift an outcome-blind SNV manifest from hg19 to hg38 with UCSC liftOver."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
    "hg19ToHg38.over.chain.gz"
)
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variants_hg19_tsv", type=Path)
    parser.add_argument("lift_over_binary", type=Path)
    parser.add_argument("chain_file", type=Path)
    parser.add_argument("variants_hg38_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    args = parser.parse_args()

    with args.variants_hg19_tsv.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    source_by_id = {row["variant_id"]: row for row in source_rows}

    with tempfile.TemporaryDirectory(prefix="alphagenome-liftover-") as temp:
        temp_path = Path(temp)
        input_bed = temp_path / "input.bed"
        output_bed = temp_path / "output.bed"
        unmapped_bed = temp_path / "unmapped.bed"
        with input_bed.open("w") as handle:
            for row in source_rows:
                start = int(row["position_1based"]) - 1
                handle.write(
                    f"{row['chromosome']}\t{start}\t{start + 1}\t"
                    f"{row['variant_id']}\t0\t+\n"
                )
        completed = subprocess.run(
            [
                str(args.lift_over_binary),
                str(input_bed),
                str(args.chain_file),
                str(output_bed),
                str(unmapped_bed),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        mappings: dict[str, list[list[str]]] = collections.defaultdict(list)
        with output_bed.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                mappings[fields[3]].append(fields)

    lifted: list[dict[str, str | int]] = []
    nonunique = 0
    for source_id, records in mappings.items():
        if len(records) != 1:
            nonunique += 1
            continue
        fields = records[0]
        source = source_by_id[source_id]
        strand = fields[5]
        reference = source["reference"]
        alternate = source["alternate"]
        if strand == "-":
            reference = reference.translate(COMPLEMENT)
            alternate = alternate.translate(COMPLEMENT)
        position = int(fields[1]) + 1
        variant_id = f"{fields[0]}:{position}:{reference}:{alternate}"
        lifted.append(
            {
                "variant_id": variant_id,
                "chromosome": fields[0],
                "position_1based": position,
                "reference": reference,
                "alternate": alternate,
                "source_variant_id_hg19": source_id,
                "source_position_1based_hg19": int(source["position_1based"]),
                "liftover_strand": strand,
                "selection_tier": source["selection_tier"],
                **{
                    key: source[key]
                    for key in (
                        "block_id_1mb",
                        "block_sha256",
                        "frozen_split",
                        "selection_rank",
                        "selection_sha256",
                        "measurement_mask",
                        "siraj_category_set",
                        "stratum_id",
                        "stratum_population",
                        "stratum_target",
                        "stratum_rank",
                    )
                    if key in source
                },
            }
        )
    lifted.sort(key=lambda row: str(row["source_variant_id_hg19"]))

    args.variants_hg38_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.variants_hg38_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(lifted[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(lifted)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_assembly": "hg19 / GRCh37",
        "target_assembly": "hg38 / GRCh38",
        "source_variants": len(source_rows),
        "uniquely_lifted_variants": len(lifted),
        "unmapped_variants": len(source_rows) - len(mappings),
        "nonunique_mappings_excluded": nonunique,
        "reverse_strand_mappings": sum(
            row["liftover_strand"] == "-" for row in lifted
        ),
        "chain_source": CHAIN_URL,
        "chain_sha256": file_sha256(args.chain_file),
        "lift_over_binary_sha256": file_sha256(args.lift_over_binary),
        "lift_over_stderr": completed.stderr.strip(),
        "mpra_functional_outcomes_accessed": False,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"source={len(source_rows)} unique_lifted={len(lifted)} "
        f"unmapped={summary['unmapped_variants']} nonunique={nonunique}"
    )


if __name__ == "__main__":
    main()
