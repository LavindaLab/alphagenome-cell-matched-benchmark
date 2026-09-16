#!/usr/bin/env python3
"""Freeze canonical SNVs in the prespecified 90% confirmation blocks.

Only variant identity and cell-measurement missingness are retained from Siraj
Table 4. Functional values and calls are never parsed or written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from inventory_pilot_strata import (
    ACTIVITY_COLUMNS,
    BLOCK_SALT,
    NS,
    column_number,
    in_gate,
    is_snv,
    sha256,
    shared_strings,
    string_cell,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table4", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    args = parser.parse_args()

    selected_columns = {1, *ACTIVITY_COLUMNS.values()}
    records = []
    with zipfile.ZipFile(args.table4) as archive:
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
                variant_cell = cells.get(1)
                variant = (
                    string_cell(variant_cell, strings)
                    if variant_cell is not None else None
                )
                if is_snv(variant) and not in_gate(variant):
                    present = []
                    for cell_name, column in ACTIVITY_COLUMNS.items():
                        cell = cells.get(column)
                        if cell is None:
                            continue
                        if (
                            cell.attrib.get("t") != "s"
                            or string_cell(cell, strings) != "NA"
                        ):
                            present.append(cell_name)
                    chromosome, position, reference, alternate = variant.split(":")
                    block_id = f"{chromosome}:{(int(position) - 1) // 1_000_000}"
                    block_hash = hashlib.sha256(
                        f"{BLOCK_SALT}\0{block_id}".encode()
                    ).hexdigest()
                    records.append({
                        "variant_id": variant,
                        "chromosome": chromosome,
                        "position_1based": int(position),
                        "reference": reference,
                        "alternate": alternate,
                        "block_id_1mb": block_id,
                        "block_sha256": block_hash,
                        "frozen_split": "confirmation",
                        "measurement_mask": "+".join(present) if present else "NONE",
                        "selection_tier": "confirmation",
                        "selection_sha256": hashlib.sha256(variant.encode()).hexdigest(),
                    })
                row.clear()

    records.sort(key=lambda item: item["variant_id"])
    if len({row["variant_id"] for row in records}) != len(records):
        raise ValueError("duplicate source variant identifiers")
    for rank, row in enumerate(records, start=1):
        row["selection_rank"] = rank
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)

    payload = "\n".join(row["variant_id"] for row in records) + "\n"
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind freeze of all canonical SNVs in confirmation blocks",
        "protocol": "PROTOCOL_v0.1.md",
        "table4_sha256": sha256(args.table4),
        "block_salt": BLOCK_SALT,
        "split_rule": "all canonical SNVs outside the lowest 10% of salted 1 Mb block-hash space",
        "confirmation_canonical_snvs": len(records),
        "unique_confirmation_blocks": len({row["block_id_1mb"] for row in records}),
        "measurement_mask_counts": dict(sorted(Counter(row["measurement_mask"] for row in records).items())),
        "source_variant_set_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "table4_fields_used": ["Variant", "cell-specific activity-estimate missingness only"],
        "numeric_mpra_values_retained": False,
        "mpra_activity_calls_accessed": False,
        "mpra_emvar_calls_accessed": False,
        "mpra_skew_values_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"confirmation_snvs={len(records)} "
        f"blocks={manifest['unique_confirmation_blocks']}"
    )


if __name__ == "__main__":
    main()
