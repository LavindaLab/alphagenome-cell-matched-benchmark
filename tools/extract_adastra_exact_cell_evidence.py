#!/usr/bin/env python3
"""Freeze exact-cell ADASTRA allele-specific binding evidence for confirmation SNVs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


CONTEXTS = {
    "A549__lung_carcinoma_": "A549",
    "HepG2__hepatoblastoma_": "HepG2",
    "K562__myelogenous_leukemia_": "K562",
    "HCT-116__colon_carcinoma_": "HCT116",
    "SK-N-SH__neuroblastoma_": "SK-N-SH",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("adastra_zip", type=Path)
    parser.add_argument("manifest_gzip", type=Path)
    parser.add_argument("output_gzip", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("--igvf-accession", default="IGVFFI5943XCOS")
    args = parser.parse_args()
    with gzip.open(args.manifest_gzip, "rt", newline="") as handle:
        manifest_ids = {row["variant_id"] for row in csv.DictReader(handle, delimiter="\t")}
    output = []
    files = 0
    archive_rows = 0
    with zipfile.ZipFile(args.adastra_zip) as archive:
        for name in archive.namelist():
            if not name.startswith("release/data/") or not name.endswith(".tsv"):
                continue
            context = next((cell for token, cell in CONTEXTS.items() if f"@{token}.tsv" in name), None)
            if context is None:
                continue
            files += 1
            tf = Path(name).name.split("@", 1)[0]
            with archive.open(name) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8")
                reader = csv.DictReader(text, delimiter="\t")
                for row in reader:
                    archive_rows += 1
                    position = int(float(row["pos"]))
                    variant = f"{row['#chr']}:{position}:{row['ref']}:{row['alt']}"
                    if variant not in manifest_ids:
                        continue
                    fdr_ref = float(row["fdrp_bh_ref"])
                    fdr_alt = float(row["fdrp_bh_alt"])
                    effect_ref = float(row["es_mean_ref"])
                    effect_alt = float(row["es_mean_alt"])
                    output.append({
                        "variant_id_hg38": variant,
                        "cell_context": context,
                        "transcription_factor": tf,
                        "binding_alt_minus_ref": effect_alt - effect_ref,
                        "fdr_ref": fdr_ref,
                        "fdr_alt": fdr_alt,
                        "asb_fdr_0_05": str(min(fdr_ref, fdr_alt) <= 0.05).lower(),
                        "total_cover": row["total_cover"],
                        "n_aggregated": row["n_aggregated"],
                        "novel_in_cltf_release": row["novel"].lower(),
                        "source_file": Path(name).name,
                    })
    output.sort(key=lambda row: (row["variant_id_hg38"], row["cell_context"], row["transcription_factor"]))
    args.output_gzip.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output_gzip, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    significant = [row for row in output if row["asb_fdr_0_05"] == "true"]
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "ADASTRA CL-TF BillCipher release (July 2022), distributed by IGVF",
        "igvf_accession": args.igvf_accession,
        "igvf_source_url": f"https://data.igvf.org/reference-files/{args.igvf_accession}/",
        "archive_sha256": sha256(args.adastra_zip),
        "manifest_sha256": sha256(args.manifest_gzip),
        "exact_cell_files": files,
        "archive_rows_scanned": archive_rows,
        "overlap_rows": len(output),
        "overlap_unique_variant_cells": len({(row["variant_id_hg38"], row["cell_context"]) for row in output}),
        "significant_fdr_0_05_rows": len(significant),
        "significant_unique_variant_cells": len({(row["variant_id_hg38"], row["cell_context"]) for row in significant}),
        "significant_rows_by_cell": dict(sorted(Counter(row["cell_context"] for row in significant).items())),
        "direction_definition": "sign(es_mean_alt - es_mean_ref); consensus across significant TFs per variant-cell is frozen for sensitivity analysis",
        "output_sha256": sha256(args.output_gzip),
        "evidence_tier": "protocol tier 2: exact-cell allele-specific TF binding",
        "confirmation_mpra_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        f"files={files} overlap_rows={len(output)} significant={len(significant)} "
        f"variant_cells={receipt['significant_unique_variant_cells']}"
    )


if __name__ == "__main__":
    main()
