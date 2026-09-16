#!/usr/bin/env python3
"""Inventory outcome-blind strata for the locked 10,000-SNV pilot.

Only the Table 4 variant identifier and missingness of the five activity-estimate
columns are retained. Numeric activity values and all activity/emVar/skew calls
are neither parsed into Python values nor written. Table 1 contributes only its
variant identifier, Category, and Variant Selection Type fields.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
BLOCK_SALT = "AGATLAS-SIRAJ-2026-V01"
CELLS = ("A549", "HEPG2", "K562", "SKNSH", "HCT116")
ACTIVITY_COLUMNS = dict(zip(CELLS, (24, 25, 26, 27, 28), strict=True))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    values: list[str] = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == f"{NS}si":
                values.append("".join(node.text or "" for node in element.iter(f"{NS}t")))
                element.clear()
    return values


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    number = 0
    for letter in letters:
        number = number * 26 + ord(letter) - 64
    return number


def string_cell(cell: ET.Element, strings: list[str]) -> str | None:
    value = cell.find(f"{NS}v")
    if value is None:
        return None
    if cell.attrib.get("t") == "s":
        return strings[int(value.text)]
    return value.text


def is_snv(variant: str | None) -> bool:
    if not variant:
        return False
    parts = variant.split(":")
    return (
        len(parts) == 4
        and len(parts[2]) == len(parts[3]) == 1
        and parts[2] in "ACGT"
        and parts[3] in "ACGT"
    )


def in_gate(variant: str) -> bool:
    chromosome, position, _, _ = variant.split(":")
    block_id = f"{chromosome}:{(int(position) - 1) // 1_000_000}"
    digest = hashlib.sha256(f"{BLOCK_SALT}\0{block_id}".encode()).digest()
    return int.from_bytes(digest) * 10 < 2**256


def table4_masks(path: Path) -> dict[str, str]:
    masks: dict[str, str] = {}
    selected = {1, *ACTIVITY_COLUMNS.values()}
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{NS}row" or int(row.attrib["r"]) < 3:
                    continue
                cells = {
                    column_number(cell.attrib["r"]): cell
                    for cell in row.findall(f"{NS}c")
                    if column_number(cell.attrib["r"]) in selected
                }
                variant_cell = cells.get(1)
                variant = string_cell(variant_cell, strings) if variant_cell is not None else None
                if is_snv(variant) and in_gate(variant):
                    present = []
                    for cell_name, column in ACTIVITY_COLUMNS.items():
                        cell = cells.get(column)
                        if cell is None:
                            continue
                        # Numeric content is classified by cell type only. For
                        # shared strings, inspect solely whether the sentinel is NA.
                        if cell.attrib.get("t") != "s" or string_cell(cell, strings) != "NA":
                            present.append(cell_name)
                    masks[variant] = "+".join(present) if present else "NONE"
                row.clear()
    return masks


def table1_metadata(path: Path, gate_variants: set[str]):
    categories: dict[str, set[str]] = collections.defaultdict(set)
    selection_types: dict[str, set[str]] = collections.defaultdict(set)
    selected = {1, 12, 17}
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{NS}row" or int(row.attrib["r"]) < 3:
                    continue
                values: dict[int, str | None] = {}
                for cell in row.findall(f"{NS}c"):
                    column = column_number(cell.attrib["r"])
                    if column in selected:
                        values[column] = string_cell(cell, strings)
                variant = values.get(1)
                if variant in gate_variants:
                    categories[variant].add(values.get(12) or "MISSING")
                    selection_types[variant].add(values.get(17) or "MISSING")
                row.clear()
    return categories, selection_types


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table4", type=Path)
    parser.add_argument("table1", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    masks = table4_masks(args.table4)
    categories, selection_types = table1_metadata(args.table1, set(masks))

    def joined(values: set[str] | None) -> str:
        return "|".join(sorted(values or {"NOT_IN_TABLE1"}))

    mask_counts = collections.Counter(masks.values())
    category_counts = collections.Counter(joined(categories.get(v)) for v in masks)
    type_counts = collections.Counter(joined(selection_types.get(v)) for v in masks)
    joint_counts = collections.Counter(
        (v.split(":", 1)[0], masks[v], joined(categories.get(v)), joined(selection_types.get(v)))
        for v in masks
    )
    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind inventory for frozen 10,000-SNV pilot stratification",
        "table4_sha256": sha256(args.table4),
        "table1_sha256": sha256(args.table1),
        "table4_fields_retained": ["Variant", "cell-specific activity-estimate missingness only"],
        "table1_fields_retained": ["Variant", "Category", "Variant Selection Type"],
        "numeric_mpra_values_retained": False,
        "mpra_activity_calls_accessed": False,
        "mpra_emvar_calls_accessed": False,
        "mpra_skew_values_accessed": False,
        "gate_canonical_snvs": len(masks),
        "measurement_mask_counts": dict(sorted(mask_counts.items())),
        "category_set_counts": dict(sorted(category_counts.items())),
        "variant_selection_type_set_counts": dict(sorted(type_counts.items())),
        "joint_strata_count": len(joint_counts),
        "joint_strata_size_distribution": dict(
            sorted(collections.Counter(joint_counts.values()).items())
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"gate_snvs={len(masks)} masks={len(mask_counts)} "
        f"categories={len(category_counts)} selection_types={len(type_counts)} "
        f"joint_strata={len(joint_counts)}"
    )


if __name__ == "__main__":
    main()
