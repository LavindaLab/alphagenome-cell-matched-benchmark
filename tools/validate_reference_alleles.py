#!/usr/bin/env python3
"""Validate 1-based SNV REF alleles against the UCSC hg38 sequence API."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


API_ROOT = "https://api.genome.ucsc.edu/getData/sequence"


def validate(
    row: dict[str, str], retries: int, assembly: str
) -> dict[str, str | int]:
    chromosome = row["chromosome"]
    position = int(row["position_1based"])
    query = urllib.parse.urlencode(
        {
            "genome": assembly,
            "chrom": chromosome,
            "start": position - 1,
            "end": position,
        }
    )
    url = f"{API_ROOT}?{query}"
    last_error = ""
    for attempt in range(1, retries + 2):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "alphagenome-metadata-feasibility/0.1"},
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = json.load(response)
            observed = str(payload.get("dna", "")).upper()
            expected = row["reference"].upper()
            status = "match" if observed == expected else "mismatch"
            return {
                **row,
                "hg38_observed_reference": observed,
                "reference_status": status,
                "attempts": attempt,
                "error": ".",
            }
        except (OSError, ValueError, urllib.error.HTTPError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt <= retries:
                time.sleep(min(2 ** (attempt - 1), 4))
    return {
        **row,
        "hg38_observed_reference": "",
        "reference_status": "retrieval_error",
        "attempts": retries + 1,
        "error": last_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variants_tsv", type=Path)
    parser.add_argument("results_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--assembly", default="hg38", choices=("hg19", "hg38"))
    args = parser.parse_args()

    with args.variants_tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        results = list(
            executor.map(
                lambda row: validate(row, args.retries, args.assembly), rows
            )
        )
    results.sort(key=lambda row: row["variant_id"])

    args.results_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.results_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(results[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    statuses = Counter(str(row["reference_status"]) for row in results)
    normalized = "\n".join(
        f"{row['variant_id']}\t{row['hg38_observed_reference']}\t{row['reference_status']}"
        for row in results
    ) + "\n"
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": args.assembly,
        "source": API_ROOT,
        "coordinate_conversion": "1-based variant position to [position-1, position) API interval",
        "variants": len(results),
        "status_counts": dict(sorted(statuses.items())),
        "valid_reference_fraction": statuses["match"] / len(results),
        "normalized_result_sha256": hashlib.sha256(normalized.encode()).hexdigest(),
        "mpra_functional_outcomes_accessed": False,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"variants={len(results)} match={statuses['match']} "
        f"mismatch={statuses['mismatch']} errors={statuses['retrieval_error']}"
    )


if __name__ == "__main__":
    main()
