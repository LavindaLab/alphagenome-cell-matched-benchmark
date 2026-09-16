#!/usr/bin/env python3
"""Merge an outcome-independent AVI recovery pass into the frozen baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def status_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=Path)
    parser.add_argument("recovery_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    base_scores = {
        row["variant_id_hg38"]: row
        for row in (json.loads(line) for line in (args.base_dir / "avi_scores.jsonl").read_text().splitlines())
    }
    recovery_scores = {
        row["variant_id_hg38"]: row
        for row in (json.loads(line) for line in (args.recovery_dir / "avi_scores.jsonl").read_text().splitlines())
    }
    base_status = {row["variant_id_hg38"]: row for row in status_rows(args.base_dir / "avi_status.tsv")}
    recovery_status = {row["variant_id_hg38"]: row for row in status_rows(args.recovery_dir / "avi_status.tsv")}
    if not set(recovery_status) <= {key for key, row in base_status.items() if row["status"] == "incomplete"}:
        raise ValueError("recovery contains IDs outside base incomplete set")
    merged_scores = {**base_scores, **recovery_scores}
    for variant, row in recovery_status.items():
        if row["status"] == "complete":
            base_status[variant] = {
                **row,
                "attempts": str(int(base_status[variant]["attempts"]) + int(row["attempts"])),
            }
        else:
            base_status[variant] = {
                **row,
                "attempts": str(int(base_status[variant]["attempts"]) + int(row["attempts"])),
            }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(row, sort_keys=True) for _, row in sorted(merged_scores.items())]
    score_bytes = ("\n".join(lines) + "\n").encode()
    (args.output_dir / "avi_scores.jsonl").write_bytes(score_bytes)
    statuses = [row for _, row in sorted(base_status.items())]
    with (args.output_dir / "avi_status.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(statuses[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(statuses)
    counts = Counter(row["status"] for row in statuses)
    base_summary = json.loads((args.base_dir / "avi_summary.json").read_text())
    recovery_summary = json.loads((args.recovery_dir / "avi_summary.json").read_text())
    merged_summary = {
        **base_summary,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(
            float(base_summary["elapsed_seconds"])
            + float(recovery_summary["elapsed_seconds"]),
            3,
        ),
        "status_counts": dict(sorted(counts.items())),
        "complete_fraction": counts["complete"] / len(statuses),
        "variant_request_attempts": sum(int(row["attempts"]) for row in statuses),
        "scores": len(merged_scores),
        "scores_sha256": hashlib.sha256(score_bytes).hexdigest(),
        "recovery_passes": 1,
        "recovery_selected_without_mpra_outcomes": True,
    }
    (args.output_dir / "avi_summary.json").write_text(
        json.dumps(merged_summary, indent=2, sort_keys=True) + "\n"
    )
    receipt = {
        "merged_utc": datetime.now(timezone.utc).isoformat(),
        "base_scores": len(base_scores),
        "recovery_requested": len(recovery_status),
        "recovery_scores": len(recovery_scores),
        "merged_scores": len(merged_scores),
        "status_counts": dict(sorted(counts.items())),
        "complete_fraction": counts["complete"] / len(statuses),
        "scores_sha256": hashlib.sha256(score_bytes).hexdigest(),
        "recovery_selected_without_mpra_outcomes": True,
        "scientific_rules_mutated": False,
    }
    (args.output_dir / "merge_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"merged_complete={counts['complete']} incomplete={counts['incomplete']}")


if __name__ == "__main__":
    main()
