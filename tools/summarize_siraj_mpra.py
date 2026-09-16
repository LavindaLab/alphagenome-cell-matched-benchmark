#!/usr/bin/env python3
"""Compute cohort-size facts from Siraj et al. supplementary XLSX tables.

The parser uses only the Python standard library and streams worksheet XML so
the 1 GB uncompressed table is never materialized in memory.
"""

from __future__ import annotations

import collections
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    values: list[str] = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == f"{NS}si":
                values.append("".join(node.text or "" for node in element.iter(f"{NS}t")))
                element.clear()
    return values


def column_number(cell_reference: str) -> int:
    letters = re.match(r"[A-Z]+", cell_reference).group(0)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - 64
    return result


def cell_value(cell: ET.Element, strings: list[str]) -> str | None:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{NS}t"))
    value_element = cell.find(f"{NS}v")
    if value_element is None:
        return None
    if cell_type == "s":
        return strings[int(value_element.text)]
    return value_element.text


def iter_selected_rows(
    archive: zipfile.ZipFile, strings: list[str], selected_columns: set[int]
):
    with archive.open("xl/worksheets/sheet1.xml") as handle:
        for _, row in ET.iterparse(handle, events=("end",)):
            if row.tag != f"{NS}row":
                continue
            row_number = int(row.attrib["r"])
            if row_number >= 3:
                values: dict[int, str | None] = {}
                for cell in row.findall(f"{NS}c"):
                    column = column_number(cell.attrib["r"])
                    if column in selected_columns:
                        values[column] = cell_value(cell, strings)
                yield values
            row.clear()


def is_snv(variant: str | None) -> bool:
    if not variant:
        return False
    parts = variant.split(":")
    return (
        len(parts) == 4
        and len(parts[2]) == 1
        and len(parts[3]) == 1
        and parts[2] in "ACGT"
        and parts[3] in "ACGT"
    )


def summarize_table4(path: Path) -> dict[str, object]:
    activity_columns = {"A549": 24, "HEPG2": 25, "K562": 26, "SKNSH": 27, "HCT116": 28}
    active_columns = {"A549": 29, "HEPG2": 30, "K562": 31, "SKNSH": 32, "HCT116": 33}
    skew_columns = {"A549": 39, "HEPG2": 40, "K562": 41, "SKNSH": 42, "HCT116": 43}
    emvar_columns = {"A549": 44, "HEPG2": 45, "K562": 46, "SKNSH": 47, "HCT116": 48}
    selected = {1, 2, 3, *activity_columns.values(), *active_columns.values(), *skew_columns.values(), *emvar_columns.values()}
    totals = collections.Counter()
    per_cell = {cell: collections.Counter() for cell in activity_columns}
    with zipfile.ZipFile(path) as archive:
        strings = load_shared_strings(archive)
        for row in iter_selected_rows(archive, strings, selected):
            variant_snv = is_snv(row.get(1))
            totals["variants"] += 1
            totals["snvs" if variant_snv else "indels_or_noncanonical"] += 1
            if row.get(2) == "1":
                totals["active_any"] += 1
                if variant_snv:
                    totals["active_any_snvs"] += 1
            if row.get(3) == "1":
                totals["emvar_any"] += 1
                if variant_snv:
                    totals["emvar_any_snvs"] += 1
            for cell in activity_columns:
                counts = per_cell[cell]
                measured = row.get(activity_columns[cell]) not in (None, "NA")
                active = row.get(active_columns[cell]) == "1"
                emvar = row.get(emvar_columns[cell]) == "1"
                skew_present = row.get(skew_columns[cell]) not in (None, "NA")
                if measured:
                    counts["measured"] += 1
                    if variant_snv:
                        counts["measured_snvs"] += 1
                if skew_present:
                    counts["skew_estimate"] += 1
                    if variant_snv:
                        counts["skew_estimate_snvs"] += 1
                if active:
                    counts["active"] += 1
                    if variant_snv:
                        counts["active_snvs"] += 1
                if emvar:
                    counts["emvar"] += 1
                    if variant_snv:
                        counts["emvar_snvs"] += 1
                if active and not emvar:
                    counts["active_non_emvar"] += 1
                    if variant_snv:
                        counts["active_non_emvar_snvs"] += 1
    return {
        "file": str(path),
        "sha256": sha256(path),
        "totals": dict(totals),
        "per_cell": {cell: dict(counts) for cell, counts in per_cell.items()},
    }


def summarize_table7(path: Path) -> dict[str, object]:
    counts_by_cell = collections.Counter()
    counts_by_cohort = collections.Counter()
    snvs_by_cell = collections.Counter()
    variants: set[str] = set()
    snvs: set[str] = set()
    with zipfile.ZipFile(path) as archive:
        strings = load_shared_strings(archive)
        for row in iter_selected_rows(archive, strings, {1, 2, 3}):
            variant = row.get(1)
            cell = row.get(2) or "missing"
            cohort = row.get(3) or "missing"
            counts_by_cell[cell] += 1
            counts_by_cohort[cohort] += 1
            if variant:
                variants.add(variant)
                if is_snv(variant):
                    snvs.add(variant)
                    snvs_by_cell[cell] += 1
    return {
        "file": str(path),
        "sha256": sha256(path),
        "rows": sum(counts_by_cell.values()),
        "unique_variants": len(variants),
        "unique_snvs": len(snvs),
        "rows_by_cell": dict(counts_by_cell),
        "snv_rows_by_cell": dict(snvs_by_cell),
        "rows_by_cohort": dict(counts_by_cohort),
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: summarize_siraj_mpra.py TABLE4.xlsx TABLE7.xlsx")
    result = {
        "table4": summarize_table4(Path(sys.argv[1])),
        "table7": summarize_table7(Path(sys.argv[2])),
    }
    json.dump(result, sys.stdout, indent=2, sort_keys=True)
    print()
