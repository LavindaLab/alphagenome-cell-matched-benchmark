#!/usr/bin/env python3
"""Select an outcome-blind, genome-spanning 100-SNV metadata sample.

Only ``xl/sharedStrings.xml`` is opened inside the MPRA workbook. Worksheet
cells, including all functional outcome columns, are never read.
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
VARIANT_RE = re.compile(
    r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X)):([1-9][0-9]*):([ACGT]):([ACGT])$"
)
CHROMOSOMES = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX"])
DEFAULT_SALT = "AGATLAS-METADATA-100-2026-V01"
PRE_OUTCOME_COMMIT = "032af87ddd8bc1489007659281bb64cf786a4d85"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selection_digest(salt: str, variant: str) -> str:
    return hashlib.sha256(f"{salt}\0{variant}".encode()).hexdigest()


def extract_variant_strings(path: Path) -> set[str]:
    variants: set[str] = set()
    with zipfile.ZipFile(path) as archive:
        if "xl/sharedStrings.xml" not in archive.namelist():
            raise ValueError("workbook has no xl/sharedStrings.xml")
        with archive.open("xl/sharedStrings.xml") as handle:
            for _, element in ET.iterparse(handle, events=("end",)):
                if element.tag != f"{NS}si":
                    continue
                value = "".join(
                    node.text or "" for node in element.iter(f"{NS}t")
                )
                match = VARIANT_RE.fullmatch(value)
                if match and match.group(3) != match.group(4):
                    variants.add(value)
                element.clear()
    return variants


def select_variants(
    variants: set[str], salt: str, per_chromosome: int = 4, total: int = 100
) -> list[dict[str, str | int]]:
    by_chromosome: dict[str, list[str]] = collections.defaultdict(list)
    for variant in variants:
        by_chromosome[variant.split(":", 1)[0]].append(variant)

    selected: list[tuple[str, str]] = []
    for chromosome in CHROMOSOMES:
        ranked = sorted(
            by_chromosome[chromosome],
            key=lambda item: (selection_digest(salt, item), item),
        )
        if len(ranked) < per_chromosome:
            raise ValueError(
                f"only {len(ranked)} canonical SNVs found on {chromosome}"
            )
        selected.extend((item, "chromosome_quota") for item in ranked[:per_chromosome])

    selected_ids = {item for item, _ in selected}
    remaining = sorted(
        variants - selected_ids,
        key=lambda item: (selection_digest(salt, item), item),
    )
    selected.extend((item, "global_fill") for item in remaining[: total - len(selected)])
    if len(selected) != total:
        raise ValueError(f"selected {len(selected)} variants; expected {total}")

    rows: list[dict[str, str | int]] = []
    for variant, tier in sorted(
        selected,
        key=lambda pair: (
            CHROMOSOMES.index(pair[0].split(":", 1)[0]),
            int(pair[0].split(":")[1]),
            pair[0],
        ),
    ):
        chromosome, position, reference, alternate = variant.split(":")
        rows.append(
            {
                "variant_id": variant,
                "chromosome": chromosome,
                "position_1based": int(position),
                "reference": reference,
                "alternate": alternate,
                "selection_tier": tier,
                "selection_sha256": selection_digest(salt, variant),
            }
        )
    return rows


def sample_sha256(rows: list[dict[str, str | int]]) -> str:
    canonical = "\n".join(sorted(str(row["variant_id"]) for row in rows)) + "\n"
    return hashlib.sha256(canonical.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_xlsx", type=Path)
    parser.add_argument("variants_tsv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    parser.add_argument("--salt", default=DEFAULT_SALT)
    args = parser.parse_args()

    variants = extract_variant_strings(args.source_xlsx)
    rows = select_variants(variants, args.salt)
    args.variants_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.variants_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    chromosome_counts = collections.Counter(str(row["chromosome"]) for row in rows)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "outcome-blind AlphaGenome Atlas metadata/API feasibility check",
        "pre_outcome_commit": PRE_OUTCOME_COMMIT,
        "source_name": args.source_xlsx.name,
        "source_sha256": file_sha256(args.source_xlsx),
        "xlsx_members_read_for_selection": ["xl/sharedStrings.xml"],
        "worksheet_cells_read_for_selection": False,
        "mpra_functional_outcomes_accessed_for_selection": False,
        "selection_salt": args.salt,
        "selection_algorithm": (
            "four lowest salted SHA-256 canonical SNVs per chr1-22/chrX, "
            "then eight globally lowest unused canonical SNVs"
        ),
        "eligible_unique_canonical_snvs": len(variants),
        "selected_variants": len(rows),
        "selected_by_chromosome": dict(sorted(chromosome_counts.items())),
        "sample_sha256": sample_sha256(rows),
    }
    args.manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"selected={len(rows)} eligible={len(variants)} "
        f"sample_sha256={manifest['sample_sha256']}"
    )


if __name__ == "__main__":
    main()
