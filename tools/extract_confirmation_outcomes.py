#!/usr/bin/env python3
"""Extract held-out MPRA outcomes only after all confirmation freezes pass."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
CELLS = {
    "A549": {"activity": 24, "active": 29, "skew": 39, "emvar": 44},
    "HepG2": {"activity": 25, "active": 30, "skew": 40, "emvar": 45},
    "K562": {"activity": 26, "active": 31, "skew": 41, "emvar": 46},
    "SK-N-SH": {"activity": 27, "active": 32, "skew": 42, "emvar": 47},
    "HCT116": {"activity": 28, "active": 33, "skew": 43, "emvar": 48},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    number = 0
    for letter in letters:
        number = number * 26 + ord(letter) - 64
    return number


def shared_strings(archive: zipfile.ZipFile):
    values = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == f"{NS}si":
                values.append("".join(node.text or "" for node in element.iter(f"{NS}t")))
                element.clear()
    return values


def value(cell, strings):
    if cell is None:
        return None
    node = cell.find(f"{NS}v")
    if node is None:
        return None
    return strings[int(node.text)] if cell.attrib.get("t") == "s" else node.text


def number(raw):
    return None if raw in (None, "", "NA") else float(raw)


def binary(raw):
    if raw in (None, "", "NA"):
        return None
    if raw not in {"0", "1"}:
        raise ValueError(f"unexpected binary value: {raw}")
    return raw == "1"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table4", type=Path)
    parser.add_argument("manifest_gzip", type=Path)
    parser.add_argument("authorization_json", type=Path)
    parser.add_argument("outcomes_gzip", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    authorization = json.loads(args.authorization_json.read_text())
    required_flags = (
        "confirmation_outcome_access_authorized",
        "prediction_retrieval_inspected",
        "analysis_implementation_committed",
        "endogenous_evidence_frozen",
        "novelty_refresh_complete",
    )
    missing = [flag for flag in required_flags if authorization.get(flag) is not True]
    if missing:
        raise ValueError(f"confirmation outcome access is not authorized: {missing}")
    with gzip.open(args.manifest_gzip, "rt", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != authorization["expected_confirmation_variants"]:
        raise ValueError("manifest count differs from authorization")
    if sha256(args.manifest_gzip) != authorization["manifest_gzip_sha256"]:
        raise ValueError("manifest checksum differs from authorization")
    final_by_source = {row["source_variant_id_hg19"]: row for row in manifest}
    if len(final_by_source) != len(manifest):
        raise ValueError("duplicate source identifiers")

    selected_columns = {1, *(field for columns in CELLS.values() for field in columns.values())}
    found = set()
    counts = collections.defaultdict(collections.Counter)
    args.outcomes_gzip.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_variant_id_hg19", "variant_id_hg38", "cell_context",
        "block_id_1mb_hg19", "measured", "activity_log2fc", "active",
        "enhancer_active", "log2skew_alt_minus_ref", "emvar",
    ]
    output_rows = 0
    with gzip.open(args.outcomes_gzip, "wt", newline="") as output, zipfile.ZipFile(args.table4) as archive:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{NS}row" or int(row.attrib["r"]) < 3:
                    continue
                cells = {
                    column_number(cell.attrib["r"]): cell
                    for cell in row.findall(f"{NS}c")
                    if column_number(cell.attrib["r"]) in selected_columns
                }
                variant = value(cells.get(1), strings)
                if variant not in final_by_source:
                    row.clear()
                    continue
                found.add(variant)
                meta = final_by_source[variant]
                for cell_name, columns in CELLS.items():
                    activity = number(value(cells.get(columns["activity"]), strings))
                    skew = number(value(cells.get(columns["skew"]), strings))
                    active = binary(value(cells.get(columns["active"]), strings))
                    emvar = binary(value(cells.get(columns["emvar"]), strings))
                    measured = activity is not None and skew is not None
                    item = {
                        "source_variant_id_hg19": variant,
                        "variant_id_hg38": meta["variant_id"],
                        "cell_context": cell_name,
                        "block_id_1mb_hg19": meta["block_id_1mb"],
                        "measured": str(measured).lower(),
                        "activity_log2fc": "" if activity is None else activity,
                        "active": "" if active is None else str(active).lower(),
                        "enhancer_active": str(bool(active and activity is not None and activity > 0)).lower(),
                        "log2skew_alt_minus_ref": "" if skew is None else skew,
                        "emvar": "NA" if emvar is None else str(emvar).lower(),
                    }
                    writer.writerow(item)
                    output_rows += 1
                    for field in ("measured", "active", "enhancer_active", "emvar"):
                        if item[field] == "true":
                            counts[cell_name][field] += 1
                row.clear()
    if set(final_by_source) != found:
        raise ValueError(f"confirmation variants absent from Table 4: {len(set(final_by_source) - found)}")
    receipt = {
        "extracted_utc": datetime.now(timezone.utc).isoformat(),
        "authorization_sha256": sha256(args.authorization_json),
        "table4_sha256": sha256(args.table4),
        "manifest_sha256": sha256(args.manifest_gzip),
        "confirmation_variants": len(manifest),
        "outcome_rows": output_rows,
        "per_cell_counts": {cell: dict(sorted(value.items())) for cell, value in sorted(counts.items())},
        "gate_outcomes_read": False,
        "confirmation_outcomes_extracted_after_freeze": True,
        "scientific_rules_mutated": False,
        "output_sha256": sha256(args.outcomes_gzip),
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"confirmation_variants={len(manifest)} outcome_rows={output_rows}")


if __name__ == "__main__":
    main()
