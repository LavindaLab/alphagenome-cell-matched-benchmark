#!/usr/bin/env python3
"""Run an outcome-blind AlphaGenome Atlas feasibility/engineering check.

This program accepts only the frozen variant manifest. It never reads or joins
the MPRA workbook or any functional outcome labels.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CONTEXTS = {
    "EFO:0001086": "A549",
    "EFO:0001187": "HepG2",
    "EFO:0002067": "K562",
    "EFO:0002824": "HCT116",
    "EFO:0003072": "SK-N-SH",
}
PRE_OUTCOME_COMMIT = "032af87ddd8bc1489007659281bb64cf786a4d85"
FORBIDDEN_INPUT_COLUMN = re.compile(
    r"emvar|skew|activity|effect|fdr|p.?value|label|outcome", re.IGNORECASE
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def load_variants(path: Path, expected_variants: int) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        forbidden = [name for name in headers if FORBIDDEN_INPUT_COLUMN.search(name)]
        if forbidden:
            raise ValueError(f"outcome-like columns are forbidden: {forbidden}")
        required = {
            "variant_id",
            "chromosome",
            "position_1based",
            "reference",
            "alternate",
        }
        if not required.issubset(headers):
            raise ValueError(f"missing columns: {sorted(required - set(headers))}")
        rows = list(reader)
    if len(rows) != expected_variants:
        raise ValueError(
            f"check requires exactly {expected_variants} variants; found {len(rows)}"
        )
    if len({row["variant_id"] for row in rows}) != len(rows):
        raise ValueError("variant manifest contains duplicate variant IDs")
    return rows


def sdk_commit(sdk_path: Path | None) -> str | None:
    if sdk_path is None:
        return None
    result = subprocess.run(
        ["git", "-C", str(sdk_path), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def import_sdk(sdk_path: Path | None):
    if sdk_path is not None:
        sys.path.insert(0, str(sdk_path / "src"))
    import alphagenome  # pylint: disable=import-outside-toplevel
    from alphagenome.atlas import atlas  # pylint: disable=import-outside-toplevel
    from alphagenome.data import genome  # pylint: disable=import-outside-toplevel
    from alphagenome.models import (  # pylint: disable=import-outside-toplevel
        variant_scorers,
    )

    return atlas, genome, variant_scorers, alphagenome.__version__


def metadata_records(frame) -> list[dict[str, Any]]:
    result = []
    for index, row in frame.iterrows():
        result.append(
            {
                "track_index": str(index),
                **{str(key): json_safe(value) for key, value in row.items()},
            }
        )
    return result


def select_required_tracks(records: list[dict[str, Any]], modality: str):
    selected = []
    for record in records:
        curie = record.get("ontology_curie")
        if curie not in CONTEXTS:
            continue
        if modality == "H3K27ac":
            mark = str(record.get("histone_mark") or "")
            name = str(record.get("name") or "")
            if mark.lower() != "h3k27ac" and "h3k27ac" not in name.lower():
                continue
        selected.append({**record, "cell_context": CONTEXTS[str(curie)]})
    return selected


def normalize_scores(
    row: dict[str, str], result: dict[str, Any], scorers: dict[str, str]
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    for modality, scorer_name in scorers.items():
        if scorer_name not in result:
            continue
        data = result[scorer_name]
        records = metadata_records(data.var)
        required = select_required_tracks(records, modality)
        quantiles = data.layers.get("quantiles")
        for record in required:
            column = list(data.var.index.map(str)).index(record["track_index"])
            raw_score = json_safe(data.X[0, column])
            quantile_score = (
                json_safe(quantiles[0, column]) if quantiles is not None else None
            )
            normalized.append(
                {
                    "variant_id_hg38": row["variant_id"],
                    "source_variant_id_hg19": row.get(
                        "source_variant_id_hg19", ""
                    ),
                    "ontology_curie": record["ontology_curie"],
                    "cell_context": record["cell_context"],
                    "modality": modality,
                    "scorer_name": scorer_name,
                    "track_index": record["track_index"],
                    "track_name": record.get("name"),
                    "assay_title": record.get("Assay title"),
                    "biosample_name": record.get("biosample_name"),
                    "histone_mark": record.get("histone_mark"),
                    "raw_score": raw_score,
                    "quantile_score": quantile_score,
                }
            )
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variants_hg38_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--sdk-path", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--expected-variants", type=int, default=100)
    parser.add_argument("--max-workers", type=int, default=5)
    parser.add_argument("--retries", type=int, default=2)
    args = parser.parse_args()

    rows = load_variants(args.variants_hg38_tsv, args.expected_variants)
    atlas, genome, variant_scorers, sdk_version = import_sdk(args.sdk_path)
    scorers = {
        "DNASE": "DNASE",
        "H3K27ac": "CHIP_HISTONE",
    }
    base_model_scorer_representations = {
        "DNASE": str(variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]),
        "H3K27ac": str(
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE"]
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    preflight = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "schema, coordinate, allele, cell-context, and retrieval validation only",
        "pre_outcome_commit": PRE_OUTCOME_COMMIT,
        "variant_manifest_sha256": file_sha256(args.variants_hg38_tsv),
        "variants": len(rows),
        "planned_variant_requests_without_retries": len(rows),
        "planned_metadata_requests": 1,
        "ontology_terms": CONTEXTS,
        "requested_scorers": scorers,
        "base_model_scorer_representations": base_model_scorer_representations,
        "sdk_version": sdk_version,
        "sdk_commit": sdk_commit(args.sdk_path),
        "mpra_workbook_read": False,
        "mpra_functional_outcomes_accessed": False,
        "outcome_aware_development": False,
        "scientific_rules_mutated": False,
    }
    (args.output_dir / "atlas_preflight.json").write_text(
        json.dumps(preflight, indent=2, sort_keys=True) + "\n"
    )
    if args.preflight_only:
        print(
            f"preflight=pass variants={len(rows)} "
            f"planned_variant_requests={len(rows)}"
        )
        return

    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise SystemExit(
            "ALPHA_GENOME_API_KEY is not set; preflight was saved and no API request was made"
        )
    client = atlas.create(api_key, timeout=30)
    scorer_metadata = client.scorer_metadata()
    scorer_inventory = {}
    for scorer_name, info in scorer_metadata.items():
        records = metadata_records(info.track_metadata)
        context_counts = Counter(
            str(record.get("ontology_curie"))
            for record in records
            if record.get("ontology_curie") in CONTEXTS
        )
        scorer_inventory[scorer_name] = {
            "is_signed": info.is_signed,
            "track_count": len(records),
            "track_metadata_columns": list(info.track_metadata.columns.map(str)),
            "required_context_track_counts": dict(sorted(context_counts.items())),
        }
    (args.output_dir / "atlas_scorer_inventory.json").write_text(
        json.dumps(scorer_inventory, indent=2, sort_keys=True) + "\n"
    )
    missing_scorers = sorted(set(scorers.values()) - set(scorer_metadata))
    if missing_scorers:
        raise RuntimeError(f"required Atlas scorers absent: {missing_scorers}")

    metadata_output: dict[str, Any] = {}
    for modality, scorer_name in scorers.items():
        info = scorer_metadata[scorer_name]
        all_records = metadata_records(info.track_metadata)
        required_records = select_required_tracks(all_records, modality)
        contexts_found = sorted(
            {str(record["ontology_curie"]) for record in required_records}
        )
        metadata_output[modality] = {
            "scorer_name": scorer_name,
            "is_signed": info.is_signed,
            "track_metadata_columns": list(info.track_metadata.columns.map(str)),
            "required_tracks": required_records,
            "contexts_found": contexts_found,
            "all_contexts_present": contexts_found == sorted(CONTEXTS),
        }
    (args.output_dir / "atlas_scorer_metadata.json").write_text(
        json.dumps(metadata_output, indent=2, sort_keys=True) + "\n"
    )
    if not all(item["all_contexts_present"] for item in metadata_output.values()):
        raise RuntimeError("not all five exact cell contexts have required tracks")

    started = time.monotonic()

    def retrieve(row: dict[str, str]):
        variant = genome.Variant(
            chromosome=row["chromosome"],
            position=int(row["position_1based"]),
            reference_bases=row["reference"],
            alternate_bases=row["alternate"],
        )
        last_error = ""
        for attempt in range(1, args.retries + 2):
            try:
                response = client.query_variant(
                    variant,
                    requested_scorers=list(scorers.values()),
                    ontology_terms=list(CONTEXTS),
                )
                return row, normalize_scores(row, response, scorers), attempt, "."
            except Exception as exc:  # bounded retries are recorded in the receipt
                last_error = " ".join(f"{type(exc).__name__}: {exc}".split())
                if attempt <= args.retries:
                    time.sleep(min(2 ** (attempt - 1), 4))
        return row, [], args.retries + 1, last_error

    retrievals = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.max_workers
    ) as executor:
        futures = [executor.submit(retrieve, row) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            retrievals.append(future.result())

    score_rows: list[dict[str, Any]] = []
    status_rows = []
    for row, scores, attempts, error in retrievals:
        score_rows.extend(scores)
        complete_pairs = {
            (score["ontology_curie"], score["modality"])
            for score in scores
            if score["raw_score"] is not None and score["quantile_score"] is not None
        }
        status_rows.append(
            {
                "variant_id_hg38": row["variant_id"],
                "source_variant_id_hg19": row.get("source_variant_id_hg19", ""),
                "status": "complete" if len(complete_pairs) == 10 else "incomplete",
                "required_context_modality_pairs_present": len(complete_pairs),
                "score_rows": len(scores),
                "attempts": attempts,
                "error": error,
            }
        )
    score_rows.sort(
        key=lambda item: (
            item["variant_id_hg38"],
            item["ontology_curie"],
            item["modality"],
            item["track_index"],
        )
    )
    status_rows.sort(key=lambda item: item["variant_id_hg38"])

    normalized_lines = [json.dumps(row, sort_keys=True) for row in score_rows]
    (args.output_dir / "atlas_scores.jsonl").write_text(
        "\n".join(normalized_lines) + ("\n" if normalized_lines else "")
    )
    with (args.output_dir / "atlas_variant_status.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(status_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(status_rows)

    statuses = Counter(str(row["status"]) for row in status_rows)
    summary = {
        **preflight,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "variant_status_counts": dict(sorted(statuses.items())),
        "complete_variant_fraction": statuses["complete"] / len(status_rows),
        "variant_request_attempts": sum(int(row["attempts"]) for row in status_rows),
        "normalized_score_rows": len(score_rows),
        "normalized_scores_sha256": hashlib.sha256(
            (("\n".join(normalized_lines) + "\n") if normalized_lines else "").encode()
        ).hexdigest(),
        "api_key_persisted": False,
        "mpra_functional_outcomes_accessed": False,
        "scientific_rules_mutated": False,
    }
    (args.output_dir / "atlas_retrieval_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"variants={len(rows)} complete={statuses['complete']} "
        f"incomplete={statuses['incomplete']} requests={summary['variant_request_attempts']}"
    )


if __name__ == "__main__":
    main()
