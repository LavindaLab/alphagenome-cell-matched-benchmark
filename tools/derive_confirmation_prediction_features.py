#!/usr/bin/env python3
"""Stream confirmation cache chunks into frozen features and AVI baselines."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")


def sign(value: float) -> int:
    return (value > 0) - (value < 0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gzip_writer(path: Path):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    return raw, compressed, text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("atlas_dir", type=Path)
    parser.add_argument("features_gzip", type=Path)
    parser.add_argument("avi_gzip", type=Path)
    parser.add_argument("status_gzip", type=Path)
    parser.add_argument("summary_json", type=Path)
    args = parser.parse_args()
    retrieval = json.loads((args.atlas_dir / "retrieval_summary.json").read_text())
    chunk_paths = sorted((args.atlas_dir / "cache_chunks").glob("chunk_*.jsonl.gz"))
    if len(chunk_paths) != retrieval["chunks"]:
        raise ValueError("cache chunk count differs from retrieval summary")

    feature_fields = [
        "variant_id_hg38", "source_variant_id_hg19", "cell_context",
        "magnitude_M", "magnitude_winning_modality", "magnitude_predicted_direction",
        "concordance_eligible", "concordance_C", "concordance_predicted_direction",
        "dnase_h3k27ac_sign_disagreement", "absolute_quantile_difference",
        "dnase_exact_vs_other_median_direction_disagreement",
        "h3k27ac_exact_vs_other_median_direction_disagreement", "tie_break_sha256",
    ]
    avi_fields = ["variant_id_hg38", "source_variant_id_hg19", "scorer_name", "raw_score", "quantile_score"]
    status_fields = ["variant_id_hg38", "source_variant_id_hg19", "status", "attempts", "error"]
    handles = []
    for path in (args.features_gzip, args.avi_gzip, args.status_gzip):
        path.parent.mkdir(parents=True, exist_ok=True)
        handles.append(gzip_writer(path))
    feature_writer = csv.DictWriter(handles[0][2], fieldnames=feature_fields, delimiter="\t", lineterminator="\n")
    avi_writer = csv.DictWriter(handles[1][2], fieldnames=avi_fields, delimiter="\t", lineterminator="\n")
    status_writer = csv.DictWriter(handles[2][2], fieldnames=status_fields, delimiter="\t", lineterminator="\n")
    feature_writer.writeheader()
    avi_writer.writeheader()
    status_writer.writeheader()
    counts = collections.Counter()
    cells = collections.Counter()
    try:
        for chunk_path in chunk_paths:
            with gzip.open(chunk_path, "rt") as chunk:
                for line in chunk:
                    record = json.loads(line)
                    counts[record["status"]] += 1
                    status_writer.writerow({field: record[field] for field in status_fields})
                    if record["status"] != "complete":
                        continue
                    avi_writer.writerow({
                        "variant_id_hg38": record["variant_id_hg38"],
                        "source_variant_id_hg19": record["source_variant_id_hg19"],
                        **record["avi"],
                    })
                    grouped = collections.defaultdict(dict)
                    other = collections.defaultdict(dict)
                    for score in record["regulatory_scores"]:
                        grouped[score["cell_context"]][score["modality"]] = score
                        other[score["modality"]][score["cell_context"]] = float(score["raw_score"])
                    if set(grouped) != set(CELLS):
                        raise ValueError(f"incomplete cells for {record['variant_id_hg38']}")
                    for cell in CELLS:
                        modalities = grouped[cell]
                        if set(modalities) != {"DNASE", "H3K27ac"}:
                            raise ValueError(f"incomplete modalities for {record['variant_id_hg38']} {cell}")
                        dnase, histone = modalities["DNASE"], modalities["H3K27ac"]
                        raw_d, raw_h = float(dnase["raw_score"]), float(histone["raw_score"])
                        q_d, q_h = float(dnase["quantile_score"]), float(histone["quantile_score"])
                        winner = "DNASE" if abs(q_d) >= abs(q_h) else "H3K27ac"
                        winner_raw = raw_d if winner == "DNASE" else raw_h
                        sign_d, sign_h = sign(raw_d), sign(raw_h)
                        concordant = sign_d != 0 and sign_d == sign_h
                        other_d = sign(statistics.median(value for name, value in other["DNASE"].items() if name != cell))
                        other_h = sign(statistics.median(value for name, value in other["H3K27ac"].items() if name != cell))
                        feature_writer.writerow({
                            "variant_id_hg38": record["variant_id_hg38"],
                            "source_variant_id_hg19": record["source_variant_id_hg19"],
                            "cell_context": cell,
                            "magnitude_M": max(abs(q_d), abs(q_h)),
                            "magnitude_winning_modality": winner,
                            "magnitude_predicted_direction": sign(winner_raw),
                            "concordance_eligible": str(concordant).lower(),
                            "concordance_C": min(abs(q_d), abs(q_h)),
                            "concordance_predicted_direction": sign_d if concordant else 0,
                            "dnase_h3k27ac_sign_disagreement": str(sign_d != sign_h).lower(),
                            "absolute_quantile_difference": abs(q_d - q_h),
                            "dnase_exact_vs_other_median_direction_disagreement": str(sign_d != other_d).lower(),
                            "h3k27ac_exact_vs_other_median_direction_disagreement": str(sign_h != other_h).lower(),
                            "tie_break_sha256": hashlib.sha256(record["variant_id_hg38"].encode()).hexdigest(),
                        })
                        cells[cell] += 1
    finally:
        for raw, compressed, text_handle in handles:
            text_handle.flush()
            text_handle.close()
            if not compressed.closed:
                compressed.close()
            if not raw.closed:
                raw.close()
    if counts["complete"] != retrieval["complete"] or counts["incomplete"] != retrieval["incomplete"]:
        raise ValueError("derived status counts differ from retrieval summary")
    if any(cells[cell] != counts["complete"] for cell in CELLS):
        raise ValueError("complete variants did not yield five feature rows")
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "complete_variants": counts["complete"],
        "incomplete_variants": counts["incomplete"],
        "prediction_feature_rows": sum(cells.values()),
        "rows_by_cell": dict(cells),
        "avi_rows": counts["complete"],
        "features_sha256": sha256(args.features_gzip),
        "avi_sha256": sha256(args.avi_gzip),
        "status_sha256": sha256(args.status_gzip),
        "cache_tree_sha256": retrieval["cache_tree_sha256"],
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"complete={counts['complete']} features={sum(cells.values())} avi={counts['complete']}")


if __name__ == "__main__":
    main()
