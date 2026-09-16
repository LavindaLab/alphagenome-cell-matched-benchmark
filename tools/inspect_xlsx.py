#!/usr/bin/env python3
"""Stream a compact schema preview from an XLSX without loading its data table."""

from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
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


def compact_preview(path: Path, row_limit: int = 3) -> dict[str, object]:
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        rows: list[list[object]] = []
        dimension = None
        max_column = 0
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for event, element in ET.iterparse(handle, events=("start", "end")):
                if event == "start" and element.tag == f"{NS}dimension":
                    dimension = element.attrib.get("ref")
                if event != "end" or element.tag != f"{NS}row":
                    continue
                sparse: dict[int, object] = {}
                for cell in element.findall(f"{NS}c"):
                    index = column_number(cell.attrib["r"])
                    max_column = max(max_column, index)
                    cell_type = cell.attrib.get("t")
                    value_element = cell.find(f"{NS}v")
                    if cell_type == "inlineStr":
                        value = "".join(
                            node.text or "" for node in cell.iter(f"{NS}t")
                        )
                    elif value_element is None:
                        value = None
                    elif cell_type == "s":
                        value = strings[int(value_element.text)]
                    else:
                        value = value_element.text
                    sparse[index] = value
                if sparse:
                    rows.append([sparse.get(i) for i in range(1, max(sparse) + 1)])
                element.clear()
                if len(rows) >= row_limit:
                    break
        return {
            "file": str(path),
            "size_bytes": path.stat().st_size,
            "dimension": dimension,
            "observed_preview_columns": max_column,
            "rows": rows,
        }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: inspect_xlsx.py WORKBOOK [WORKBOOK ...]")
    json.dump(
        [compact_preview(Path(name)) for name in sys.argv[1:]],
        sys.stdout,
        indent=2,
    )
    print()
