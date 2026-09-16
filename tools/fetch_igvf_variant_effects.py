#!/usr/bin/env python3
"""Snapshot official IGVF Variant-EFFECTS edges before confirmation outcomes."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


BASE = "https://api.catalogkg.igvf.org/api/variants/genes"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_jsonl_gz", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("--page-size", type=int, default=500)
    parser.add_argument("--max-pages", type=int, default=1000)
    args = parser.parse_args()
    args.output_jsonl_gz.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    pages = 0
    started = time.monotonic()
    with gzip.open(args.output_jsonl_gz, "wt") as output:
        for page in range(args.max_pages):
            query = urllib.parse.urlencode({
                "method": "Variant-EFFECTS",
                "source": "IGVF",
                "verbose": "true",
                "limit": args.page_size,
                "page": page,
            })
            request = urllib.request.Request(f"{BASE}?{query}", headers={"User-Agent": "AlphaGenome-external-validation/0.1"})
            with urllib.request.urlopen(request, timeout=120) as response:
                records = json.load(response)
            if not isinstance(records, list):
                raise ValueError("IGVF endpoint did not return a list")
            if not records:
                break
            for record in records:
                output.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
            total += len(records)
            pages += 1
            if len(records) < args.page_size:
                break
        else:
            raise RuntimeError("max-pages reached before an empty/partial page")
    receipt = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source": "official IGVF Catalog production API",
        "endpoint": BASE,
        "query": {"method": "Variant-EFFECTS", "source": "IGVF", "verbose": True},
        "pages": pages,
        "records": total,
        "page_size": args.page_size,
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "output_sha256": sha256(args.output_jsonl_gz),
        "confirmation_mpra_outcomes_accessed": False,
        "selection_tier": "protocol tier 1: endogenous genome editing with allele-level molecular outcome",
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"pages={pages} records={total}")


if __name__ == "__main__":
    main()
