#!/usr/bin/env python3
"""Validate SNV REF alleles locally against an official UCSC twoBit assembly."""

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


HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_fasta(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    name: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    result[name] = "".join(sequence).upper()
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if name is not None:
        result[name] = "".join(sequence).upper()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variants_tsv", type=Path)
    parser.add_argument("two_bit", type=Path)
    parser.add_argument("two_bit_to_fa", type=Path)
    parser.add_argument("results_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("--assembly", default="hg38")
    args = parser.parse_args()

    with args.variants_tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with tempfile.TemporaryDirectory(prefix="alphagenome-twobit-") as temp:
        temp_path = Path(temp)
        bed = temp_path / "variants.bed"
        fasta = temp_path / "variants.fa"
        with bed.open("w") as handle:
            for index, row in enumerate(rows):
                start = int(row["position_1based"]) - 1
                handle.write(f"{row['chromosome']}\t{start}\t{start + 1}\tv{index}\n")
        subprocess.run(
            [str(args.two_bit_to_fa), f"-bed={bed}", str(args.two_bit), str(fasta)],
            check=True,
            capture_output=True,
            text=True,
        )
        observed = parse_fasta(fasta)

    results = []
    for index, row in enumerate(rows):
        base = observed.get(f"v{index}", "")
        status = "match" if base == row["reference"].upper() else "mismatch"
        results.append(
            {
                **row,
                f"{args.assembly}_observed_reference": base,
                "reference_status": status,
                "attempts": 1,
                "error": ".",
            }
        )
    results.sort(key=lambda row: row["variant_id"])
    args.results_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.results_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    counts = collections.Counter(row["reference_status"] for row in results)
    normalized = "\n".join(
        f"{row['variant_id']}\t{row[f'{args.assembly}_observed_reference']}\t{row['reference_status']}"
        for row in results
    ) + "\n"
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": args.assembly,
        "source": HG38_URL,
        "two_bit_sha256": sha256(args.two_bit),
        "two_bit_to_fa_sha256": sha256(args.two_bit_to_fa),
        "coordinate_conversion": "1-based variant position to [position-1, position) BED interval",
        "variants": len(results),
        "status_counts": dict(sorted(counts.items())),
        "valid_reference_fraction": counts["match"] / len(results),
        "normalized_result_sha256": hashlib.sha256(normalized.encode()).hexdigest(),
        "failed_sandbox_api_receipt_preserved": True,
        "mpra_functional_outcomes_accessed": False,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"variants={len(results)} match={counts['match']} mismatch={counts['mismatch']}")


if __name__ == "__main__":
    main()
