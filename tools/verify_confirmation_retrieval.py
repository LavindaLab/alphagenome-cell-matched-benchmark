#!/usr/bin/env python3
"""Independently verify frozen confirmation cache and derived predictions."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
MODALITIES = ("DNASE", "H3K27ac")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_bytes(records) -> bytes:
    payload = "\n".join(json.dumps(row, sort_keys=True, separators=(",", ":")) for row in records)
    return (payload + ("\n" if payload else "")).encode()


def sign(value: float) -> int:
    return (value > 0) - (value < 0)


def tsv_count_and_ids(path: Path, id_field: str):
    digest = hashlib.sha256()
    count = 0
    identifiers = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            count += 1
            identifiers.add(row[id_field])
            digest.update((row[id_field] + "\n").encode())
    return count, identifiers, digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest_gzip", type=Path)
    parser.add_argument("atlas_dir", type=Path)
    parser.add_argument("features_gzip", type=Path)
    parser.add_argument("avi_gzip", type=Path)
    parser.add_argument("status_gzip", type=Path)
    parser.add_argument("derivation_summary", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()

    summary = json.loads((args.atlas_dir / "retrieval_summary.json").read_text())
    metadata_path = args.atlas_dir / "atlas_scorer_metadata.json"
    metadata_sha = sha256(metadata_path)
    chunk_paths = sorted((args.atlas_dir / "cache_chunks").glob("chunk_*.jsonl.gz"))
    receipt_paths = sorted((args.atlas_dir / "cache_chunks").glob("chunk_*.receipt.json"))
    checks = []

    def check(name, passed, observed=None, expected=None):
        checks.append({"check": name, "passed": bool(passed), "observed": observed, "expected": expected})

    with gzip.open(args.manifest_gzip, "rt", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    manifest_map = {row["variant_id"]: row["source_variant_id_hg19"] for row in manifest}
    check("manifest_unique", len(manifest_map) == len(manifest), len(manifest_map), len(manifest))
    check("manifest_sha", sha256(args.manifest_gzip) == summary["manifest_sha256"], sha256(args.manifest_gzip), summary["manifest_sha256"])
    check("chunk_file_count", len(chunk_paths) == summary["chunks"], len(chunk_paths), summary["chunks"])
    check("chunk_receipt_count", len(receipt_paths) == summary["chunks"], len(receipt_paths), summary["chunks"])

    statuses = Counter()
    attempts = 0
    seen = set()
    complete_ids = set()
    duplicate_ids = 0
    invalid_records = 0
    chunk_receipts = []
    for index, (chunk_path, chunk_receipt_path) in enumerate(zip(chunk_paths, receipt_paths, strict=True), start=1):
        receipt = json.loads(chunk_receipt_path.read_text())
        chunk_receipts.append(receipt)
        compressed_sha = sha256(chunk_path)
        with gzip.open(chunk_path, "rb") as handle:
            uncompressed = handle.read()
        lines = uncompressed.splitlines()
        check(f"chunk_{index:04d}_receipt", (
            receipt["chunk"] == index
            and receipt["passed"] is True
            and receipt["records"] == len(lines)
            and receipt["file_sha256"] == compressed_sha
            and receipt["uncompressed_sha256"] == hashlib.sha256(uncompressed).hexdigest()
        ))
        for line in lines:
            record = json.loads(line)
            variant = record["variant_id_hg38"]
            if variant in seen:
                duplicate_ids += 1
            seen.add(variant)
            statuses[record["status"]] += 1
            attempts += int(record["attempts"])
            valid = (
                variant in manifest_map
                and record["source_variant_id_hg19"] == manifest_map[variant]
                and record["scorer_metadata_sha256"] == metadata_sha
                and record["mpra_functional_outcomes_accessed"] is False
            )
            if record["status"] == "complete":
                complete_ids.add(variant)
                regulatory = record["regulatory_scores"]
                pairs = {(item["cell_context"], item["modality"]) for item in regulatory}
                numeric = [item[field] for item in regulatory for field in ("raw_score", "quantile_score")]
                numeric.extend(record["avi"][field] for field in ("raw_score", "quantile_score"))
                response = {"regulatory_scores": regulatory, "avi": record["avi"]}
                valid = valid and (
                    len(regulatory) == 10
                    and pairs == {(cell, modality) for cell in CELLS for modality in MODALITIES}
                    and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in numeric)
                    and record["response_sha256"] == hashlib.sha256(canonical_bytes([response])).hexdigest()
                )
            else:
                valid = valid and record["status"] == "incomplete"
            invalid_records += not valid

    tree_payload = "\n".join(f"{item['chunk']:04d}\t{item['file_sha256']}" for item in chunk_receipts) + "\n"
    check("cache_tree_sha", hashlib.sha256(tree_payload.encode()).hexdigest() == summary["cache_tree_sha256"])
    check("cache_ids_exact", seen == set(manifest_map), len(seen), len(manifest_map))
    check("cache_no_duplicates", duplicate_ids == 0, duplicate_ids, 0)
    check("cache_record_schema_and_hashes", invalid_records == 0, invalid_records, 0)
    check("summary_complete", statuses["complete"] == summary["complete"], statuses["complete"], summary["complete"])
    check("summary_incomplete", statuses["incomplete"] == summary["incomplete"], statuses["incomplete"], summary["incomplete"])
    check("summary_attempts", attempts == summary["request_attempts"], attempts, summary["request_attempts"])
    check("retrieval_completeness_gate", summary["complete_fraction"] >= 0.99, summary["complete_fraction"], 0.99)
    check("outcome_blind_summary", summary["mpra_functional_outcomes_accessed"] is False)

    derivation = json.loads(args.derivation_summary.read_text())
    feature_count, feature_ids, _ = tsv_count_and_ids(args.features_gzip, "variant_id_hg38")
    avi_count, avi_ids, _ = tsv_count_and_ids(args.avi_gzip, "variant_id_hg38")
    status_count, status_ids, _ = tsv_count_and_ids(args.status_gzip, "variant_id_hg38")
    check("feature_row_count", feature_count == statuses["complete"] * len(CELLS), feature_count, statuses["complete"] * len(CELLS))
    check("feature_variant_ids", feature_ids == complete_ids, len(feature_ids), len(complete_ids))
    check("avi_row_count", avi_count == statuses["complete"], avi_count, statuses["complete"])
    check("avi_variant_ids", avi_ids == complete_ids, len(avi_ids), len(complete_ids))
    check("status_row_count", status_count == len(manifest), status_count, len(manifest))
    check("status_variant_ids", status_ids == set(manifest_map), len(status_ids), len(manifest_map))
    check("features_sha", sha256(args.features_gzip) == derivation["features_sha256"])
    check("avi_sha", sha256(args.avi_gzip) == derivation["avi_sha256"])
    check("status_sha", sha256(args.status_gzip) == derivation["status_sha256"])
    check("derivation_cache_tree", derivation["cache_tree_sha256"] == summary["cache_tree_sha256"])
    check("derivation_outcome_blind", derivation["mpra_functional_outcomes_accessed"] is False)

    derived_mismatches = 0
    derived_records = 0
    with (
        gzip.open(args.features_gzip, "rt", newline="") as feature_handle,
        gzip.open(args.avi_gzip, "rt", newline="") as avi_handle,
        gzip.open(args.status_gzip, "rt", newline="") as status_handle,
    ):
        feature_reader = iter(csv.DictReader(feature_handle, delimiter="\t"))
        avi_reader = iter(csv.DictReader(avi_handle, delimiter="\t"))
        status_reader = iter(csv.DictReader(status_handle, delimiter="\t"))
        for chunk_path in chunk_paths:
            with gzip.open(chunk_path, "rt") as chunk_handle:
                for line in chunk_handle:
                    record = json.loads(line)
                    status_row = next(status_reader)
                    derived_records += 1
                    status_expected = {
                        "variant_id_hg38": record["variant_id_hg38"],
                        "source_variant_id_hg19": record["source_variant_id_hg19"],
                        "status": record["status"],
                        "attempts": str(record["attempts"]),
                        "error": record["error"],
                    }
                    derived_mismatches += status_row != status_expected
                    if record["status"] != "complete":
                        continue
                    avi_row = next(avi_reader)
                    derived_mismatches += not (
                        avi_row["variant_id_hg38"] == record["variant_id_hg38"]
                        and avi_row["source_variant_id_hg19"] == record["source_variant_id_hg19"]
                        and avi_row["scorer_name"] == "AVI_SCORE"
                        and float(avi_row["raw_score"]) == float(record["avi"]["raw_score"])
                        and float(avi_row["quantile_score"]) == float(record["avi"]["quantile_score"])
                    )
                    grouped = {(item["cell_context"], item["modality"]): item for item in record["regulatory_scores"]}
                    for cell in CELLS:
                        feature = next(feature_reader)
                        dnase = grouped[(cell, "DNASE")]
                        histone = grouped[(cell, "H3K27ac")]
                        raw_d, raw_h = float(dnase["raw_score"]), float(histone["raw_score"])
                        q_d, q_h = float(dnase["quantile_score"]), float(histone["quantile_score"])
                        winner = "DNASE" if abs(q_d) >= abs(q_h) else "H3K27ac"
                        sign_d, sign_h = sign(raw_d), sign(raw_h)
                        other_d = sign(statistics.median(float(grouped[(other, "DNASE")]["raw_score"]) for other in CELLS if other != cell))
                        other_h = sign(statistics.median(float(grouped[(other, "H3K27ac")]["raw_score"]) for other in CELLS if other != cell))
                        expected = {
                            "variant_id_hg38": record["variant_id_hg38"],
                            "source_variant_id_hg19": record["source_variant_id_hg19"],
                            "cell_context": cell,
                            "magnitude_M": max(abs(q_d), abs(q_h)),
                            "magnitude_winning_modality": winner,
                            "magnitude_predicted_direction": sign(raw_d if winner == "DNASE" else raw_h),
                            "concordance_eligible": str(sign_d != 0 and sign_d == sign_h).lower(),
                            "concordance_C": min(abs(q_d), abs(q_h)),
                            "concordance_predicted_direction": sign_d if sign_d != 0 and sign_d == sign_h else 0,
                            "dnase_h3k27ac_sign_disagreement": str(sign_d != sign_h).lower(),
                            "absolute_quantile_difference": abs(q_d - q_h),
                            "dnase_exact_vs_other_median_direction_disagreement": str(sign_d != other_d).lower(),
                            "h3k27ac_exact_vs_other_median_direction_disagreement": str(sign_h != other_h).lower(),
                            "tie_break_sha256": hashlib.sha256(record["variant_id_hg38"].encode()).hexdigest(),
                        }
                        for key, value in expected.items():
                            observed = feature[key]
                            if isinstance(value, float):
                                derived_mismatches += float(observed) != value
                            else:
                                derived_mismatches += str(observed) != str(value)
        derived_mismatches += sum(1 for _ in feature_reader)
        derived_mismatches += sum(1 for _ in avi_reader)
        derived_mismatches += sum(1 for _ in status_reader)
    check("derived_values_match_cache", derived_mismatches == 0, derived_mismatches, 0)
    check("derived_records_match_manifest", derived_records == len(manifest), derived_records, len(manifest))

    passed = all(item["passed"] for item in checks)
    receipt = {
        "inspected_utc": datetime.now(timezone.utc).isoformat(),
        "passed": passed,
        "checks_passed": sum(item["passed"] for item in checks),
        "checks_total": len(checks),
        "checks": checks,
        "manifest_variants": len(manifest),
        "complete": statuses["complete"],
        "incomplete": statuses["incomplete"],
        "complete_fraction": statuses["complete"] / len(manifest),
        "feature_rows": feature_count,
        "cache_tree_sha256": summary["cache_tree_sha256"],
        "features_sha256": sha256(args.features_gzip),
        "avi_sha256": sha256(args.avi_gzip),
        "status_sha256": sha256(args.status_gzip),
        "confirmation_mpra_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    args.receipt_json.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"passed={passed} checks={receipt['checks_passed']}/{receipt['checks_total']} complete={statuses['complete']}/{len(manifest)}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
