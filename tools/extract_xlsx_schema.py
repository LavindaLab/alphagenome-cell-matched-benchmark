#!/usr/bin/env python3
"""Extract only title/header rows from an XLSX worksheet.

The parser stops before the first data row, ensuring that no row-level outcome
values are accessed during schema inspection.
"""

from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - 64
    return result


def load_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    values: list[str] = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == f"{NS}si":
                values.append(
                    "".join(node.text or "" for node in element.iter(f"{NS}t"))
                )
                element.clear()
    return values


def cell_value(cell: ET.Element, strings: list[str]) -> str | None:
    value = cell.find(f"{NS}v")
    if value is None:
        return None
    if cell.attrib.get("t") == "s":
        return strings[int(value.text)]
    return value.text


def extract(path: Path) -> dict[str, object]:
    with zipfile.ZipFile(path) as archive:
        strings = load_shared_strings(archive)
        title = None
        headers: dict[int, str | None] = {}
        dimension = None
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for event, element in ET.iterparse(handle, events=("start", "end")):
                if event == "start" and element.tag == f"{NS}dimension":
                    dimension = element.attrib.get("ref")
                if event != "end" or element.tag != f"{NS}row":
                    continue
                row_number = int(element.attrib["r"])
                if row_number > 2:
                    break
                cells = element.findall(f"{NS}c")
                if row_number == 1 and cells:
                    title = cell_value(cells[0], strings)
                elif row_number == 2:
                    for cell in cells:
                        headers[column_number(cell.attrib["r"])] = cell_value(
                            cell, strings
                        )
                element.clear()
    return {
        "source_name": path.name,
        "dimension": dimension,
        "title": title,
        "headers": [headers[index] for index in range(1, max(headers) + 1)],
        "rows_read": [1, 2],
        "data_rows_read": False,
        "mpra_functional_outcomes_accessed": False,
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: extract_xlsx_schema.py WORKBOOK OUTPUT.json")
    result = extract(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + "\n")
    print(f"columns={len(result['headers'])} data_rows_read=false")
