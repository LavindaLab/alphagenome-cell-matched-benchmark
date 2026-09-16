#!/usr/bin/env python3
"""Extract MPRA outcomes only for the frozen gate-pilot variants.

The caller must run this only after committing the Atlas prediction freeze.
Non-selected worksheet rows are identified by Variant and discarded before any
functional cell is converted to a Python value.
"""

from __future__ import annotations

import argparse
import collections
import csv
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


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    values = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == f"{NS}si":
                values.append("".join(node.text or "" for node in element.iter(f"{NS}t")))
                element.clear()
    return values


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - 64
    return result


def value(cell: ET.Element | None, strings: list[str]) -> str | None:
    if cell is None:
        return None
    node = cell.find(f"{NS}v")
    if node is None:
        return None
    if cell.attrib.get("t") == "s":
        return strings[int(node.text)]
    return node.text


def number(raw: str | None) -> float | None:
    if raw in (None, "", "NA"):
        return None
    return float(raw)


def binary(raw: str | None) -> bool | None:
    if raw in (None, "", "NA"):
        return None
    if raw not in {"0", "1"}:
        raise ValueError(f"unexpected binary value: {raw}")
    return raw == "1"


def extract_table4(path: Path, final_by_source: dict[str, dict[str, str]]):
    output = []
    found = set()
    selected_columns = {1, *(field for columns in CELLS.values() for field in columns.values())}
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{NS}row" or int(row.attrib["r"]) < 3:
                    continue
                cells = {}
                variant_cell = None
                for cell in row.findall(f"{NS}c"):
                    column = column_number(cell.attrib["r"])
                    if column == 1:
                        variant_cell = cell
                    if column in selected_columns:
                        cells[column] = cell
                variant = value(variant_cell, strings)
                if variant not in final_by_source:
                    row.clear()
                    continue
                found.add(variant)
                metadata = final_by_source[variant]
                for cell_name, columns in CELLS.items():
                    activity = number(value(cells.get(columns["activity"]), strings))
                    skew = number(value(cells.get(columns["skew"]), strings))
                    active = binary(value(cells.get(columns["active"]), strings))
                    emvar = binary(value(cells.get(columns["emvar"]), strings))
                    measured = activity is not None and skew is not None
                    output.append(
                        {
                            "source_variant_id_hg19": variant,
                            "variant_id_hg38": metadata["variant_id"],
                            "cell_context": cell_name,
                            "block_id_1mb_hg19": metadata["block_id_1mb"],
                            "stratum_id": metadata["stratum_id"],
                            "measured": str(measured).lower(),
                            "activity_log2fc": "" if activity is None else activity,
                            "active": "" if active is None else str(active).lower(),
                            "enhancer_active": str(bool(active and activity is not None and activity > 0)).lower(),
                            "log2skew_alt_minus_ref": "" if skew is None else skew,
                            "emvar": "NA" if emvar is None else str(emvar).lower(),
                        }
                    )
                row.clear()
    missing = set(final_by_source) - found
    if missing:
        raise ValueError(f"frozen pilot variants absent from Table 4: {len(missing)}")
    return output


def extract_table1(path: Path, final_by_source: dict[str, dict[str, str]]):
    aggregated = {
        variant: {"pips": [], "categories": set(), "selection_types": set(), "cs_ids": set()}
        for variant in final_by_source
    }
    selected = {1, 10, 12, 17, 18}
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{NS}row" or int(row.attrib["r"]) < 3:
                    continue
                cells = {}
                variant_cell = None
                for cell in row.findall(f"{NS}c"):
                    column = column_number(cell.attrib["r"])
                    if column == 1:
                        variant_cell = cell
                    if column in selected:
                        cells[column] = cell
                variant = value(variant_cell, strings)
                if variant not in aggregated:
                    row.clear()
                    continue
                item = aggregated[variant]
                pip = number(value(cells.get(10), strings))
                if pip is not None:
                    item["pips"].append(pip)
                item["categories"].add(value(cells.get(12), strings) or "MISSING")
                item["selection_types"].add(value(cells.get(17), strings) or "MISSING")
                item["cs_ids"].add(value(cells.get(18), strings) or "MISSING")
                row.clear()
    rows = []
    for source_id, item in sorted(aggregated.items()):
        categories = sorted(item["categories"] or {"NOT_IN_TABLE1"})
        cohort_categories = [category for category in categories if category != "control"]
        rows.append(
            {
                "source_variant_id_hg19": source_id,
                "variant_id_hg38": final_by_source[source_id]["variant_id"],
                "category_set": "|".join(categories),
                "fine_mapped_candidate": str(bool(cohort_categories)).lower(),
                "max_pip": "" if not item["pips"] else max(item["pips"]),
                "pip_available": str(bool(item["pips"])).lower(),
                "variant_selection_type_set": "|".join(sorted(item["selection_types"] or {"NOT_IN_TABLE1"})),
                "credible_set_id_count": len(item["cs_ids"] - {"MISSING"}),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table4", type=Path)
    parser.add_argument("table1", type=Path)
    parser.add_argument("final_hg38_manifest", type=Path)
    parser.add_argument("outcomes_tsv", type=Path)
    parser.add_argument("variant_metadata_tsv", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()
    with args.final_hg38_manifest.open(newline="") as handle:
        final = list(csv.DictReader(handle, delimiter="\t"))
    if len(final) != 10_000 or any(row["frozen_split"] != "gate" for row in final):
        raise ValueError("outcome extraction is restricted to exactly 10,000 frozen gate variants")
    final_by_source = {row["source_variant_id_hg19"]: row for row in final}
    if len(final_by_source) != 10_000:
        raise ValueError("duplicate source identifiers in final manifest")
    outcomes = extract_table4(args.table4, final_by_source)
    metadata = extract_table1(args.table1, final_by_source)
    args.outcomes_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outcomes_tsv, outcomes)
    write_tsv(args.variant_metadata_tsv, metadata)
    per_cell = collections.defaultdict(collections.Counter)
    for row in outcomes:
        counts = per_cell[row["cell_context"]]
        for field in ("measured", "active", "enhancer_active", "emvar"):
            if row[field] == "true":
                counts[field] += 1
    receipt = {
        "extracted_utc": datetime.now(timezone.utc).isoformat(),
        "authorization_boundary": "gate-pilot outcomes only after committed Atlas prediction freeze 00e55d9",
        "table4_sha256": sha256(args.table4),
        "table1_sha256": sha256(args.table1),
        "final_manifest_sha256": sha256(args.final_hg38_manifest),
        "frozen_gate_variants": len(final),
        "outcome_rows": len(outcomes),
        "variant_metadata_rows": len(metadata),
        "per_cell_counts": {cell: dict(sorted(counts.items())) for cell, counts in sorted(per_cell.items())},
        "table4_fields_accessed_for_selected_gate_rows": [
            "Variant", "cell-specific activity log2FC", "cell-specific active call",
            "cell-specific log2Skew", "cell-specific emVar call",
        ],
        "confirmation_outcomes_retained_or_written": False,
        "confirmation_split_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"gate_variants={len(final)} outcome_rows={len(outcomes)} metadata_rows={len(metadata)}")


if __name__ == "__main__":
    main()
