#!/usr/bin/env python3
"""Run the frozen held-out primary, utility, discordance, and endogenous analyses."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import os
import random
import sqlite3
import statistics
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


CELLS = ("A549", "HepG2", "K562", "HCT116", "SK-N-SH")
CONFIRMATION_CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")
COVERAGES = (0.01, 0.05, 0.10, 0.25, 0.50)
BUDGETS = (96, 384, 1000, 5000)
POLICIES = ("random", "pip", "avi", "magnitude", "concordance", "concordance_block_cap")
BOOTSTRAPS = 2000
SEED = 20260913


@dataclass(slots=True)
class Obs:
    variant: str
    source: str
    chromosome: str
    position: int
    block_1mb: str
    outcome_direction: int
    abs_skew: float
    emvar: bool
    enhancer_active: bool
    active: bool
    magnitude: float
    magnitude_direction: int
    concordance_eligible: bool
    concordance: float
    concordance_direction: int
    sign_disagreement: bool
    quantile_difference: float
    dnase_other_disagreement: bool
    histone_other_disagreement: bool
    tie_hash: str
    fine_mapped: bool
    pip: float | None
    avi: float
    promoter: bool


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def truth(value: str) -> int:
    return int(value == "true")


def sign(value: float) -> int:
    return (value > 0) - (value < 0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path):
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", newline="")
    else:
        handle = path.open(newline="")
    with handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, records) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_tsv_gzip(path: Path, records) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def percentile(values, probability):
    return float(np.quantile(np.asarray(values, dtype=float), probability, method="linear"))


def batch_insert(connection, sql, iterable, batch_size=10000):
    batch = []
    count = 0
    for item in iterable:
        batch.append(item)
        if len(batch) >= batch_size:
            connection.executemany(sql, batch)
            count += len(batch)
            batch.clear()
    if batch:
        connection.executemany(sql, batch)
        count += len(batch)
    connection.commit()
    return count


def build_database(connection, args):
    connection.executescript("""
        PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE; PRAGMA cache_size=-500000;
        CREATE TABLE metadata(variant TEXT PRIMARY KEY, fine INTEGER, pip REAL);
        CREATE TABLE regions(variant TEXT PRIMARY KEY, promoter INTEGER);
        CREATE TABLE avi(variant TEXT PRIMARY KEY, score REAL);
        CREATE TABLE features(
          variant TEXT, cell TEXT, magnitude REAL, magnitude_direction INTEGER,
          concordance_eligible INTEGER, concordance REAL, concordance_direction INTEGER,
          sign_disagreement INTEGER, quantile_difference REAL,
          dnase_other_disagreement INTEGER, histone_other_disagreement INTEGER,
          tie_hash TEXT, PRIMARY KEY(variant, cell));
        CREATE TABLE outcomes(
          variant TEXT, source TEXT, cell TEXT, block_1mb TEXT, skew REAL,
          emvar INTEGER, enhancer_active INTEGER, active INTEGER,
          PRIMARY KEY(variant, cell));
    """)
    counts = {}
    counts["metadata"] = batch_insert(
        connection, "INSERT INTO metadata VALUES (?,?,?)",
        ((r["variant_id_hg38"], truth(r["fine_mapped_candidate"]), None if r["max_pip"] == "" else float(r["max_pip"])) for r in rows(args.metadata)),
    )
    counts["regions"] = batch_insert(
        connection, "INSERT INTO regions VALUES (?,?)",
        ((r["variant_id_hg38"], truth(r["promoter_proximal_2kb"])) for r in rows(args.regions)),
    )
    counts["avi"] = batch_insert(
        connection, "INSERT INTO avi VALUES (?,?)",
        ((r["variant_id_hg38"], float(r["quantile_score"])) for r in rows(args.avi)),
    )
    counts["features"] = batch_insert(
        connection, "INSERT INTO features VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        ((
            r["variant_id_hg38"], r["cell_context"], float(r["magnitude_M"]), int(r["magnitude_predicted_direction"]),
            truth(r["concordance_eligible"]), float(r["concordance_C"]), int(r["concordance_predicted_direction"]),
            truth(r["dnase_h3k27ac_sign_disagreement"]), float(r["absolute_quantile_difference"]),
            truth(r["dnase_exact_vs_other_median_direction_disagreement"]), truth(r["h3k27ac_exact_vs_other_median_direction_disagreement"]),
            r["tie_break_sha256"],
        ) for r in rows(args.features)),
    )
    counts["outcomes"] = batch_insert(
        connection, "INSERT INTO outcomes VALUES (?,?,?,?,?,?,?,?)",
        ((
            r["variant_id_hg38"], r["source_variant_id_hg19"], r["cell_context"], r["block_id_1mb_hg19"],
            float(r["log2skew_alt_minus_ref"]), truth(r["emvar"]), truth(r["enhancer_active"]), truth(r["active"]),
        ) for r in rows(args.outcomes) if r["measured"] == "true"),
    )
    connection.executescript("CREATE INDEX outcomes_cell ON outcomes(cell); CREATE INDEX features_cell ON features(cell);")
    return counts


def load_cell(connection, cell):
    query = """
      SELECT o.variant,o.source,o.block_1mb,o.skew,o.emvar,o.enhancer_active,o.active,
             f.magnitude,f.magnitude_direction,f.concordance_eligible,f.concordance,f.concordance_direction,
             f.sign_disagreement,f.quantile_difference,f.dnase_other_disagreement,f.histone_other_disagreement,f.tie_hash,
             m.fine,m.pip,a.score,r.promoter
      FROM outcomes o JOIN features f USING(variant,cell) JOIN metadata m USING(variant)
      JOIN avi a USING(variant) JOIN regions r USING(variant) WHERE o.cell=?
    """
    result = []
    for value in connection.execute(query, (cell,)):
        chromosome, position, _, _ = value[1].split(":")
        result.append(Obs(
            variant=value[0], source=value[1], chromosome=chromosome, position=int(position), block_1mb=value[2],
            outcome_direction=sign(float(value[3])), abs_skew=abs(float(value[3])), emvar=bool(value[4]),
            enhancer_active=bool(value[5]), active=bool(value[6]), magnitude=float(value[7]), magnitude_direction=int(value[8]),
            concordance_eligible=bool(value[9]), concordance=float(value[10]), concordance_direction=int(value[11]),
            sign_disagreement=bool(value[12]), quantile_difference=float(value[13]),
            dnase_other_disagreement=bool(value[14]), histone_other_disagreement=bool(value[15]), tie_hash=value[16],
            fine_mapped=bool(value[17]), pip=None if value[18] is None else float(value[18]), avi=float(value[19]), promoter=bool(value[20]),
        ))
    return result


def primary_estimate(population, coverage):
    requested = math.ceil(len(population) * coverage)
    magnitude = sorted(population, key=lambda r: (-r.magnitude, r.tie_hash))[:requested]
    concordance = sorted((r for r in population if r.concordance_eligible), key=lambda r: (-r.concordance, r.tie_hash))[:requested]
    magnitude_error = statistics.fmean(r.magnitude_direction != r.outcome_direction for r in magnitude)
    concordance_error = statistics.fmean(r.concordance_direction != r.outcome_direction for r in concordance)
    return {
        "population": len(population), "requested": requested, "magnitude_selected": len(magnitude),
        "concordance_selected": len(concordance), "magnitude_error": magnitude_error,
        "concordance_error": concordance_error, "difference": concordance_error - magnitude_error,
    }


def random_stratified(population, budget, cell):
    strata = collections.defaultdict(list)
    for item in population:
        strata[(item.chromosome, item.pip is not None)].append(item)
    requested = min(budget, len(population))
    quotas = {key: requested * len(value) / len(population) for key, value in strata.items()}
    allocation = {key: math.floor(value) for key, value in quotas.items()}
    remaining = requested - sum(allocation.values())
    order = sorted(strata, key=lambda key: (-(quotas[key] - allocation[key]), str(key)))
    for key in order[:remaining]:
        allocation[key] += 1
    selected = []
    for key, values in strata.items():
        selected.extend(sorted(values, key=lambda r: hashlib.sha256(f"AGATLAS-UTILITY-RANDOM-2026-V01\0{r.source}\0{cell}".encode()).hexdigest())[:allocation[key]])
    return selected


def select_policy(population, policy, budget, cell):
    requested = min(budget, len(population))
    if policy == "random":
        return random_stratified(population, requested, cell)
    if policy == "pip":
        return sorted((r for r in population if r.pip is not None), key=lambda r: (-r.pip, r.tie_hash))[:requested]
    if policy == "avi":
        return sorted(population, key=lambda r: (-r.avi, r.tie_hash))[:requested]
    if policy == "magnitude":
        return sorted(population, key=lambda r: (-r.magnitude, r.tie_hash))[:requested]
    concordance = sorted((r for r in population if r.concordance_eligible), key=lambda r: (-r.concordance, r.tie_hash))
    if policy == "concordance":
        return concordance[:requested]
    selected = []
    block_counts = collections.Counter()
    for item in concordance:
        if block_counts[item.block_1mb] >= 2:
            continue
        selected.append(item)
        block_counts[item.block_1mb] += 1
        if len(selected) == requested:
            break
    return selected


def block_key(item, size):
    return item.chromosome, (item.position - 1) // size


def weighted_error(errors, block_indices, weights, requested):
    row_weights = weights[block_indices]
    available = int(row_weights.sum())
    selected = min(requested, available)
    if selected == 0:
        return math.nan
    cumulative = np.cumsum(row_weights)
    end = int(np.searchsorted(cumulative, selected, side="left"))
    before = int(cumulative[end - 1]) if end else 0
    error_count = float(np.dot(errors[:end], row_weights[:end]))
    error_count += float(errors[end]) * (selected - before)
    return error_count / selected


def bootstrap_primary(primary_by_cell, size, seed):
    all_blocks = sorted({block_key(item, size) for values in primary_by_cell.values() for item in values})
    block_index = {block: index for index, block in enumerate(all_blocks)}
    by_chromosome = collections.defaultdict(list)
    for index, (chromosome, _) in enumerate(all_blocks):
        by_chromosome[chromosome].append(index)
    prepared = {}
    for cell, population in primary_by_cell.items():
        magnitude = sorted(population, key=lambda r: (-r.magnitude, r.tie_hash))
        concordance = sorted((r for r in population if r.concordance_eligible), key=lambda r: (-r.concordance, r.tie_hash))
        prepared[cell] = {
            "population_blocks": np.asarray([block_index[block_key(r, size)] for r in population], dtype=np.int32),
            "magnitude_blocks": np.asarray([block_index[block_key(r, size)] for r in magnitude], dtype=np.int32),
            "magnitude_errors": np.asarray([r.magnitude_direction != r.outcome_direction for r in magnitude], dtype=np.int8),
            "concordance_blocks": np.asarray([block_index[block_key(r, size)] for r in concordance], dtype=np.int32),
            "concordance_errors": np.asarray([r.concordance_direction != r.outcome_direction for r in concordance], dtype=np.int8),
        }
    rng = random.Random(seed)
    draws = []
    for replicate in range(BOOTSTRAPS):
        weights = np.zeros(len(all_blocks), dtype=np.int16)
        for indices in by_chromosome.values():
            for _ in indices:
                weights[rng.choice(indices)] += 1
        differences = {}
        for cell, data in prepared.items():
            sample_n = int(weights[data["population_blocks"]].sum())
            requested = math.ceil(sample_n * 0.10)
            mag = weighted_error(data["magnitude_errors"], data["magnitude_blocks"], weights, requested)
            con = weighted_error(data["concordance_errors"], data["concordance_blocks"], weights, requested)
            differences[cell] = con - mag
        draws.append({
            "replicate": replicate + 1, "block_size": size,
            "confirmation_macro_difference": statistics.fmean(differences[cell] for cell in CONFIRMATION_CELLS),
            "five_cell_descriptive_macro_difference": statistics.fmean(differences.values()),
            **{f"{cell}_difference": differences[cell] for cell in CELLS},
        })
    return draws


def adastra_sensitivity(connection, path):
    grouped = collections.defaultdict(list)
    for row in rows(path):
        if row["asb_fdr_0_05"] == "true":
            grouped[(row["variant_id_hg38"], row["cell_context"])].append(sign(float(row["binding_alt_minus_ref"])))
    output = []
    for cell in CELLS:
        consensus = []
        ambiguous = 0
        for (variant, context), directions in grouped.items():
            if context != cell:
                continue
            unique = {value for value in directions if value != 0}
            if len(unique) != 1:
                ambiguous += 1
                continue
            feature = connection.execute(
                "SELECT magnitude_direction,concordance_eligible,concordance_direction FROM features WHERE variant=? AND cell=?",
                (variant, cell),
            ).fetchone()
            if feature is not None:
                consensus.append((next(iter(unique)), int(feature[0]), bool(feature[1]), int(feature[2])))
        concordant = [row for row in consensus if row[2]]
        output.append({
            "cell_context": cell,
            "consensus_variant_cells": len(consensus),
            "ambiguous_multi_tf_variant_cells": ambiguous,
            "magnitude_direction_accuracy": statistics.fmean(row[0] == row[1] for row in consensus) if consensus else math.nan,
            "concordance_eligible_variant_cells": len(concordant),
            "concordance_direction_accuracy": statistics.fmean(row[0] == row[3] for row in concordant) if concordant else math.nan,
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("outcomes", type=Path)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("features", type=Path)
    parser.add_argument("avi", type=Path)
    parser.add_argument("regions", type=Path)
    parser.add_argument("adastra", type=Path)
    parser.add_argument("authorization", type=Path)
    parser.add_argument("outcome_receipt", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    authorization = json.loads(args.authorization.read_text())
    required_flags = (
        "confirmation_outcome_access_authorized",
        "prediction_retrieval_inspected",
        "analysis_implementation_committed",
        "endogenous_evidence_frozen",
        "novelty_refresh_complete",
    )
    missing = [flag for flag in required_flags if authorization.get(flag) is not True]
    if missing:
        raise ValueError(f"confirmation analysis is not authorized: {missing}")
    expected_hashes = authorization.get("analysis_input_sha256", {})
    for name in ("metadata", "features", "avi", "regions", "adastra"):
        observed = sha256(getattr(args, name))
        if expected_hashes.get(name) != observed:
            raise ValueError(f"{name} checksum differs from authorization")
    if authorization.get("analysis_script_sha256") != sha256(Path(__file__)):
        raise ValueError("analysis script checksum differs from authorization")
    outcome_receipt = json.loads(args.outcome_receipt.read_text())
    if outcome_receipt.get("output_sha256") != sha256(args.outcomes):
        raise ValueError("outcome checksum differs from extraction receipt")
    if outcome_receipt.get("authorization_sha256") != sha256(args.authorization):
        raise ValueError("outcome receipt was not produced under this authorization")
    if outcome_receipt.get("gate_outcomes_read") is not False:
        raise ValueError("outcome receipt does not preserve the gate/confirmation firewall")
    if outcome_receipt.get("confirmation_outcomes_extracted_after_freeze") is not True:
        raise ValueError("outcomes were not extracted after the confirmation freeze")
    expected_variants = int(authorization["expected_confirmation_variants"])
    if outcome_receipt.get("confirmation_variants") != expected_variants:
        raise ValueError("outcome receipt variant count differs from authorization")
    if outcome_receipt.get("outcome_rows") != expected_variants * len(CELLS):
        raise ValueError("outcome receipt does not contain five rows per confirmation variant")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory is not empty; preserve the first frozen run")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    descriptor, database_name = tempfile.mkstemp(prefix="alphagenome-confirmation-", suffix=".sqlite")
    os.close(descriptor)
    connection = sqlite3.connect(database_name)
    try:
        import_counts = build_database(connection, args)
        if import_counts["metadata"] != expected_variants or import_counts["regions"] != expected_variants:
            raise ValueError("outcome-blind variant inputs are incomplete")
        if import_counts["avi"] < math.ceil(expected_variants * 0.99):
            raise ValueError("AVI completeness is below the frozen 99% gate")
        if import_counts["features"] != import_counts["avi"] * len(CELLS):
            raise ValueError("feature rows are not complete across five cells")
        primary_by_cell = {}
        primary_rows = []
        risk_rows = []
        utility_rows = []
        discordance_rows = []
        for cell in CELLS:
            cell_rows = load_cell(connection, cell)
            primary = [row for row in cell_rows if row.enhancer_active and row.emvar]
            primary_by_cell[cell] = primary
            for coverage in COVERAGES:
                estimate = primary_estimate(primary, coverage)
                risk_rows.append({"cell_context": cell, "coverage": coverage, **estimate})
                if coverage == 0.10:
                    primary_rows.append({
                        "cell_context": cell,
                        "confirmation_role": "confirmatory" if cell in CONFIRMATION_CELLS else "descriptive_count_rule_downgrade",
                        **estimate,
                    })
            universe = [row for row in cell_rows if row.fine_mapped]
            total_emvars = sum(row.emvar for row in universe)
            for budget in BUDGETS:
                for policy in POLICIES:
                    selected = select_policy(universe, policy, budget, cell)
                    hits = [row for row in selected if row.emvar]
                    utility_rows.append({
                        "cell_context": cell, "budget": budget, "policy": policy,
                        "eligible_population": len(universe), "selected": len(selected),
                        "emvar_hits": len(hits), "precision": len(hits) / len(selected) if selected else math.nan,
                        "cumulative_absolute_mpra_skew": sum(row.abs_skew for row in selected),
                        "independent_blocks_with_hit": len({row.block_1mb for row in hits}),
                        "promoter_proximal_selected": sum(row.promoter for row in selected),
                        "promoter_proximal_hits": sum(row.promoter for row in hits),
                        "missed_emvars": total_emvars - len(hits),
                    })
            for disagrees in (False, True):
                subset = [row for row in cell_rows if row.sign_disagreement == disagrees]
                enhancer_emvars = [row for row in subset if row.emvar and row.enhancer_active]
                discordance_rows.append({
                    "cell_context": cell,
                    "dnase_h3k27ac_sign_disagreement": str(disagrees).lower(),
                    "measured_noncoding": len(subset), "emvars": sum(row.emvar for row in subset),
                    "emvar_prevalence": statistics.fmean(row.emvar for row in subset) if subset else math.nan,
                    "enhancer_active_emvars": len(enhancer_emvars),
                    "magnitude_sign_error_among_enhancer_active_emvars": statistics.fmean(row.magnitude_direction != row.outcome_direction for row in enhancer_emvars) if enhancer_emvars else math.nan,
                })
            del cell_rows

        primary_macro = statistics.fmean(row["difference"] for row in primary_rows if row["cell_context"] in CONFIRMATION_CELLS)
        five_cell_macro = statistics.fmean(row["difference"] for row in primary_rows)
        draws_1mb = bootstrap_primary(primary_by_cell, 1_000_000, SEED)
        draws_5mb = bootstrap_primary(primary_by_cell, 5_000_000, SEED + 5)
        draws = draws_1mb + draws_5mb
        ci_1mb = [percentile([row["confirmation_macro_difference"] for row in draws_1mb], 0.025), percentile([row["confirmation_macro_difference"] for row in draws_1mb], 0.975)]
        ci_5mb = [percentile([row["confirmation_macro_difference"] for row in draws_5mb], 0.025), percentile([row["confirmation_macro_difference"] for row in draws_5mb], 0.975)]
        adastra_rows = adastra_sensitivity(connection, args.adastra)

        utility_lookup = {(row["cell_context"], int(row["budget"]), row["policy"]): row for row in utility_rows}
        utility_screens = {}
        for budget in (384, 1000):
            mag = [utility_lookup[(cell, budget, "magnitude")]["precision"] for cell in CONFIRMATION_CELLS]
            avi = [utility_lookup[(cell, budget, "avi")]["precision"] for cell in CONFIRMATION_CELLS]
            macro_mag, macro_avi = statistics.fmean(mag), statistics.fmean(avi)
            utility_screens[str(budget)] = {
                "magnitude_macro_precision": macro_mag,
                "avi_macro_precision": macro_avi,
                "absolute_difference": macro_mag - macro_avi,
                "relative_gain": macro_mag / macro_avi - 1 if macro_avi else math.inf,
                "cells_favoring_magnitude": sum(left > right for left, right in zip(mag, avi, strict=True)),
            }
        screen_a = all(
            value["absolute_difference"] >= 0.03 and value["relative_gain"] >= 0.20 and value["cells_favoring_magnitude"] >= 3
            for value in utility_screens.values()
        )
        primary_differences = [row["difference"] for row in primary_rows if row["cell_context"] in CONFIRMATION_CELLS]
        same_direction = sum((value > 0) == (primary_macro > 0) and value != 0 for value in primary_differences)
        screen_b = abs(primary_macro) >= 0.02 and (ci_1mb[0] > 0 or ci_1mb[1] < 0) and same_direction >= 3
        by_cell_discordance = collections.defaultdict(dict)
        for row in discordance_rows:
            by_cell_discordance[row["cell_context"]][row["dnase_h3k27ac_sign_disagreement"]] = row
        agreement_rr = {}
        for cell in CONFIRMATION_CELLS:
            numerator = by_cell_discordance[cell]["false"]["emvar_prevalence"]
            denominator = by_cell_discordance[cell]["true"]["emvar_prevalence"]
            agreement_rr[cell] = numerator / denominator if denominator > 0 else None
        screen_c = sum(value is not None and value >= 1.5 for value in agreement_rr.values()) >= 3 and ci_1mb[0] > -0.02
        result = {
            "generated_utc": now(),
            "primary_estimand": "four-cell macro mean sign-error(concordant)-sign-error(magnitude) at 10% coverage among noncoding enhancer-active emVars",
            "confirmation_cells": list(CONFIRMATION_CELLS),
            "a549_role": "descriptive_count_rule_downgrade",
            "primary_macro_difference": primary_macro,
            "primary_1mb_bootstrap_95_ci": ci_1mb,
            "primary_5mb_bootstrap_95_ci": ci_5mb,
            "five_cell_descriptive_macro_difference": five_cell_macro,
            "primary_counts_by_cell": {row["cell_context"]: row["population"] for row in primary_rows},
            "bootstrap_replicates": BOOTSTRAPS,
            "bootstrap_seed": SEED,
            "utility_screens": utility_screens,
            "agreement_emvar_prevalence_ratio_by_confirmation_cell": agreement_rr,
            "quantitative_value_screens": {"A_experimental_utility": screen_a, "B_concordance_policy": screen_b, "C_marker_policy_distinction": screen_c},
            "qualitative_value_gate_pending": True,
            "import_counts": import_counts,
            "scientific_rules_mutated": False,
        }
        write_tsv(args.output_dir / "primary_by_cell.tsv", primary_rows)
        write_tsv(args.output_dir / "risk_coverage.tsv", risk_rows)
        write_tsv(args.output_dir / "budget_utility.tsv", utility_rows)
        write_tsv(args.output_dir / "discordance.tsv", discordance_rows)
        write_tsv(args.output_dir / "endogenous_adastra.tsv", adastra_rows)
        write_tsv_gzip(args.output_dir / "bootstrap_draws.tsv.gz", draws)
        (args.output_dir / "confirmation_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        file_hashes = {
            path.name: sha256(path) for path in sorted(args.output_dir.iterdir()) if path.is_file()
        }
        manifest = {
            "generated_utc": now(),
            "input_sha256": {
                **{name: sha256(getattr(args, name)) for name in ("outcomes", "outcome_receipt", "metadata", "features", "avi", "regions", "adastra", "authorization")},
                "analysis_script": sha256(Path(__file__)),
            },
            "output_sha256": file_hashes,
            "confirmation_outcomes_used": True,
            "gate_outcomes_used": False,
            "scientific_rules_mutated": False,
        }
        (args.output_dir / "result_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(f"primary_difference={primary_macro:.6f} ci95=[{ci_1mb[0]:.6f},{ci_1mb[1]:.6f}] screens={result['quantitative_value_screens']}")
    finally:
        connection.close()
        Path(database_name).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
