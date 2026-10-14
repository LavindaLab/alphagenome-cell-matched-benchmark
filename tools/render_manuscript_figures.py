#!/usr/bin/env python3
"""Render manuscript figures and traceable source-data tables from committed outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Line, PolyLine, Rect, String
from reportlab.lib import colors
from reportlab.pdfgen import canvas


WIDTH = 540
HEIGHT = 610
NAVY = colors.HexColor("#24476B")
ORANGE = colors.HexColor("#D9792B")
TEAL = colors.HexColor("#2A9D8F")
PURPLE = colors.HexColor("#7251A4")
GRAY = colors.HexColor("#6B7280")
LIGHT = colors.HexColor("#E5E7EB")
DARK = colors.HexColor("#1F2937")
RED = colors.HexColor("#B24A4A")
PALETTE = {
    "magnitude": NAVY,
    "avi": ORANGE,
    "pip": PURPLE,
    "random": GRAY,
    "concordance": TEAL,
    "concordance_block_cap": RED,
}
LABELS = {
    "magnitude": "Cell-matched",
    "avi": "AVI",
    "pip": "PIP",
    "random": "Random",
    "concordance": "Concordance",
    "concordance_block_cap": "Concordance + block cap",
}
CELLS = ("HepG2", "K562", "HCT116", "SK-N-SH")
BUDGETS = (96, 384, 1000, 5000)


def read_tsv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, rows):
    fieldnames = []
    for row in rows:
        for name in row:
            if name not in fieldnames:
                fieldnames.append(name)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def text(drawing, x, y, value, size=8, color=DARK, anchor="start", bold=False):
    drawing.add(String(x, y, str(value), fontName="Helvetica-Bold" if bold else "Helvetica", fontSize=size, fillColor=color, textAnchor=anchor))


def panel(drawing, label, title, x, y):
    text(drawing, x, y, label, 11, DARK, bold=True)
    text(drawing, x + 18, y, title, 9, DARK, bold=True)


def axes(drawing, left, bottom, width, height, x_ticks, x_labels, y_min, y_max, y_ticks, y_label=None):
    drawing.add(Line(left, bottom, left, bottom + height, strokeColor=DARK, strokeWidth=0.7))
    drawing.add(Line(left, bottom, left + width, bottom, strokeColor=DARK, strokeWidth=0.7))
    for value in y_ticks:
        y = bottom + (value - y_min) / (y_max - y_min) * height
        drawing.add(Line(left, y, left + width, y, strokeColor=LIGHT, strokeWidth=0.4))
        text(drawing, left - 5, y - 2.5, f"{value:g}", 6.5, GRAY, "end")
    for position, label in zip(x_ticks, x_labels, strict=True):
        drawing.add(Line(position, bottom, position, bottom - 3, strokeColor=DARK, strokeWidth=0.5))
        text(drawing, position, bottom - 12, label, 6.5, GRAY, "middle")
    if y_label:
        text(drawing, left - 34, bottom + height / 2, y_label, 7, DARK, "middle")


def forest(drawing, rows, left, bottom, width, height, x_min, x_max, value_key, low_key, high_key, labels, zero=0.0, color=NAVY):
    drawing.add(Line(left, bottom, left + width, bottom, strokeColor=DARK, strokeWidth=0.7))
    if x_min <= zero <= x_max:
        x_zero = left + (zero - x_min) / (x_max - x_min) * width
        drawing.add(Line(x_zero, bottom, x_zero, bottom + height, strokeColor=GRAY, strokeWidth=0.7, strokeDashArray=[2, 2]))
    for tick in np.linspace(x_min, x_max, 5):
        x = left + (tick - x_min) / (x_max - x_min) * width
        drawing.add(Line(x, bottom, x, bottom - 3, strokeColor=DARK, strokeWidth=0.4))
        text(drawing, x, bottom - 12, f"{tick:.2f}", 6.2, GRAY, "middle")
    gap = height / max(len(rows), 1)
    for index, (row, label) in enumerate(zip(rows, labels, strict=True)):
        y = bottom + height - gap * (index + 0.5)
        value = float(row[value_key]); low = float(row[low_key]); high = float(row[high_key])
        x_value = left + (value - x_min) / (x_max - x_min) * width
        x_low = left + (low - x_min) / (x_max - x_min) * width
        x_high = left + (high - x_min) / (x_max - x_min) * width
        drawing.add(Line(x_low, y, x_high, y, strokeColor=color, strokeWidth=1.2))
        drawing.add(Line(x_low, y - 2, x_low, y + 2, strokeColor=color, strokeWidth=0.8))
        drawing.add(Line(x_high, y - 2, x_high, y + 2, strokeColor=color, strokeWidth=0.8))
        drawing.add(Circle(x_value, y, 2.6, fillColor=color, strokeColor=color))
        text(drawing, left - 6, y - 2.5, label, 6.7, DARK, "end")


def figure_1(project):
    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Independent, outcome-blind evaluation of AlphaGenome Atlas", 14, DARK, bold=True)
    panel(drawing, "A", "Frozen evaluation workflow", 25, 548)
    boxes = [
        ("Siraj MPRA", "304,278 variants"),
        ("Coordinate/QC", "246,674 held-out SNVs"),
        ("Atlas retrieval", "5 exact cell contexts"),
        ("Prediction freeze", "outcomes inaccessible"),
        ("Held-out analysis", "rules unchanged"),
    ]
    x_positions = [25, 130, 235, 340, 445]
    for index, ((title, subtitle), x) in enumerate(zip(boxes, x_positions, strict=True)):
        drawing.add(Rect(x, 475, 84, 45, rx=4, ry=4, fillColor=colors.white, strokeColor=NAVY, strokeWidth=1))
        text(drawing, x + 42, 502, title, 7.2, NAVY, "middle", True)
        text(drawing, x + 42, 486, subtitle, 6.2, GRAY, "middle")
        if index < len(boxes) - 1:
            drawing.add(Line(x + 84, 497, x + 101, 497, strokeColor=GRAY, strokeWidth=1))
            drawing.add(Line(x + 97, 500, x + 101, 497, strokeColor=GRAY, strokeWidth=1))
            drawing.add(Line(x + 97, 494, x + 101, 497, strokeColor=GRAY, strokeWidth=1))

    panel(drawing, "B", "Exact-context evaluation and frozen roles", 25, 435)
    cells = [
        ("A549", "EFO:0001086", "descriptive", "2,249"),
        ("HepG2", "EFO:0001187", "confirmatory", "13,045"),
        ("K562", "EFO:0002067", "confirmatory", "15,809"),
        ("HCT116", "EFO:0002824", "confirmatory", "5,866"),
        ("SK-N-SH", "EFO:0003072", "confirmatory", "13,762"),
    ]
    headers = ["MPRA cell", "Atlas ontology", "role", "directional n"]
    columns = [35, 145, 285, 430]
    for x, header in zip(columns, headers, strict=True):
        text(drawing, x, 405, header, 7, GRAY, bold=True)
    for index, values in enumerate(cells):
        y = 382 - index * 28
        if index % 2 == 0:
            drawing.add(Rect(25, y - 8, 490, 22, fillColor=colors.HexColor("#F7F8FA"), strokeColor=None))
        for x, value in zip(columns, values, strict=True):
            text(drawing, x, y, value, 7.5, DARK)

    panel(drawing, "C", "Integrity and independence checks", 25, 225)
    checks = [
        ("Atlas completeness", "246,674 / 246,674", TEAL),
        ("Prediction inspection", "273 / 273 checks", TEAL),
        ("Frozen-result inspection", "25 / 25 checks", TEAL),
        ("Manuscript analysis inspection", "100 / 100 checks", TEAL),
        ("Gate outcomes in confirmation", "0", NAVY),
        ("Duplicate variant-cell rows", "0", NAVY),
    ]
    for index, (label, value, color) in enumerate(checks):
        row = index // 2; col = index % 2
        x = 30 + col * 255; y = 183 - row * 52
        drawing.add(Rect(x, y, 230, 36, rx=3, ry=3, fillColor=colors.HexColor("#F7F8FA"), strokeColor=LIGHT))
        text(drawing, x + 10, y + 22, label, 7, GRAY)
        text(drawing, x + 220, y + 12, value, 9, color, "end", True)
    text(drawing, 25, 18, "All policy definitions, thresholds and outcomes were frozen before held-out outcome access.", 7.2, GRAY)
    return drawing


def figure_2(project, source):
    utility = read_tsv(project / "confirmation/results/budget_utility.tsv")
    contrasts = read_tsv(project / "manuscript/analysis/budget_precision_contrasts.tsv")
    controls = project / "manuscript/analysis/context-controls"
    context_rows = [
        row for row in read_tsv(controls / "wrong_cell_contrasts.tsv")
        if int(row["block_size_bp"]) == 1_000_000
    ]
    track_rows = [
        row for row in read_tsv(controls / "single_track_estimates.tsv")
        if int(row["block_size_bp"]) == 1_000_000 and row["cell_context"] == "confirmation_macro"
    ]
    macro_rows = []
    for budget in BUDGETS:
        for policy in ("random", "pip", "avi", "magnitude", "concordance"):
            values = [float(row["precision"]) for row in utility if row["cell_context"] in CELLS and int(row["budget"]) == budget and row["policy"] == policy]
            macro_rows.append({"budget": budget, "policy": policy, "macro_precision": np.mean(values)})
    write_tsv(source / "figure_2a_macro_precision.tsv", macro_rows)
    write_tsv(source / "figure_2b_context_macro.tsv", [row for row in context_rows if row["cell_context"] == "confirmation_macro"])
    write_tsv(source / "figure_2c_context_cell_contrasts.tsv", context_rows)
    write_tsv(source / "figure_2e_single_track_macro.tsv", track_rows)
    yield_rows = []
    for budget in (384, 1000):
        for policy in ("magnitude", "avi"):
            selected = [row for row in utility if row["cell_context"] in CELLS and int(row["budget"]) == budget and row["policy"] == policy]
            yield_rows.append({
                "budget": budget, "policy": policy,
                "emvar_hits": sum(int(row["emvar_hits"]) for row in selected),
                "hit_containing_1mb_blocks": sum(int(row["independent_blocks_with_hit"]) for row in selected),
            })
    write_tsv(source / "figure_2d_yield.tsv", yield_rows)

    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Cell context improves AlphaGenome regulatory-variant prioritization", 13.5, DARK, bold=True)
    panel(drawing, "A", "Macro emVar precision across fixed assay budgets", 25, 548)
    left, bottom, width, height = 65, 382, 430, 135
    x_positions = [left + (math.log10(budget) - math.log10(96)) / (math.log10(5000) - math.log10(96)) * width for budget in BUDGETS]
    axes(drawing, left, bottom, width, height, x_positions, [str(x) for x in BUDGETS], 0, 0.85, [0, 0.2, 0.4, 0.6, 0.8], "Precision")
    for policy in ("magnitude", "concordance", "avi", "pip", "random"):
        values = [next(float(row["macro_precision"]) for row in macro_rows if int(row["budget"]) == budget and row["policy"] == policy) for budget in BUDGETS]
        points = []
        for x, value in zip(x_positions, values, strict=True):
            y = bottom + value / 0.85 * height
            points.extend([x, y]); drawing.add(Circle(x, y, 2.3, fillColor=PALETTE[policy], strokeColor=PALETTE[policy]))
        drawing.add(PolyLine(points, strokeColor=PALETTE[policy], strokeWidth=1.6, fillColor=None))
    for index, policy in enumerate(("magnitude", "concordance", "avi", "pip", "random")):
        x = 70 + index * 88
        drawing.add(Line(x, 352, x + 13, 352, strokeColor=PALETTE[policy], strokeWidth=2))
        text(drawing, x + 17, 349, LABELS[policy], 6.2, DARK)
    text(drawing, 280, 365, "Assay slots per cell (log scale)", 7, DARK, "middle")

    panel(drawing, "B", "Exact-cell versus nonmatching-context mean", 25, 326)
    macro_context = [row for row in context_rows if row["cell_context"] == "confirmation_macro"]
    b_left, b_bottom, b_width, b_height = 55, 170, 195, 125
    axes(drawing, b_left, b_bottom, b_width, b_height, [], [], 0, 0.75, [0, 0.25, 0.50, 0.75], "Precision")
    for index, budget in enumerate((384, 1000)):
        row = next(item for item in macro_context if int(item["budget"]) == budget)
        center = b_left + b_width * (index + 0.5) / 2
        for offset, key, color in ((-10, "exact_cell_precision", NAVY), (10, "all_wrong_cell_mean_precision", GRAY)):
            value = float(row[key])
            drawing.add(Rect(center + offset - 7, b_bottom, 14, value / 0.75 * b_height, fillColor=color, strokeColor=None))
        text(drawing, center, b_bottom - 13, f"{budget:,}", 6.5, DARK, "middle")
    drawing.add(Rect(80, 146, 10, 7, fillColor=NAVY, strokeColor=None)); text(drawing, 95, 146, "Exact cell", 6.4)
    drawing.add(Rect(150, 146, 10, 7, fillColor=GRAY, strokeColor=None)); text(drawing, 165, 146, "Nonmatching mean", 6.4)

    panel(drawing, "C", "Exact-cell benefit by outcome cell", 285, 326)
    c_left, c_bottom, c_width, c_height = 365, 170, 145, 125
    drawing.add(Line(c_left, c_bottom, c_left + c_width, c_bottom, strokeColor=DARK, strokeWidth=0.7))
    x_zero = c_left
    drawing.add(Line(x_zero, c_bottom, x_zero, c_bottom + c_height, strokeColor=GRAY, strokeWidth=0.7, strokeDashArray=[2, 2]))
    for tick in (0, 0.1, 0.2, 0.3, 0.4):
        x = c_left + tick / 0.42 * c_width
        drawing.add(Line(x, c_bottom, x, c_bottom - 3, strokeColor=DARK, strokeWidth=0.4))
        text(drawing, x, c_bottom - 12, f"{tick:.1f}", 6.1, GRAY, "middle")
    gap = c_height / len(CELLS)
    for index, cell in enumerate(CELLS):
        y_center = c_bottom + c_height - gap * (index + 0.5)
        text(drawing, c_left - 7, y_center - 2, cell, 6.4, DARK, "end")
        for y_offset, budget, color in ((4, 384, NAVY), (-4, 1000, TEAL)):
            row = next(item for item in context_rows if item["cell_context"] == cell and int(item["budget"]) == budget)
            conv = lambda value: c_left + float(value) / 0.42 * c_width
            y = y_center + y_offset
            drawing.add(Line(conv(row["ci95_low"]), y, conv(row["ci95_high"]), y, strokeColor=color, strokeWidth=1.1))
            drawing.add(Circle(conv(row["exact_minus_all_wrong_mean"]), y, 2.2, fillColor=color, strokeColor=color))
    drawing.add(Line(373, 146, 386, 146, strokeColor=NAVY, strokeWidth=2)); text(drawing, 391, 143, "384", 6.3)
    drawing.add(Line(435, 146, 448, 146, strokeColor=TEAL, strokeWidth=2)); text(drawing, 453, 143, "1,000", 6.3)

    panel(drawing, "D", "Additional confirmed variants and genomic blocks", 25, 116)
    text(drawing, 45, 88, "384 slots", 7, GRAY, bold=True)
    text(drawing, 45, 70, "+417 emVars", 9, NAVY, bold=True)
    text(drawing, 45, 53, "+367 hit-containing 1 Mb blocks", 6.8, DARK)
    text(drawing, 165, 88, "1,000 slots", 7, GRAY, bold=True)
    text(drawing, 165, 70, "+973 emVars", 9, NAVY, bold=True)
    text(drawing, 165, 53, "+674 hit-containing 1 Mb blocks", 6.8, DARK)

    panel(drawing, "E", "Post hoc single-track macro precision", 285, 116)
    e_left, e_bottom, e_width, e_height = 325, 42, 180, 46
    track_order = (("avi", ORANGE), ("dnase_only", PURPLE), ("h3k27ac_only", TEAL), ("combined_magnitude", NAVY))
    for budget_index, budget in enumerate((384, 1000)):
        group_y = e_bottom + e_height - budget_index * 30
        text(drawing, e_left - 8, group_y - 3, f"{budget:,}", 6.2, DARK, "end")
        for track_index, (score, color) in enumerate(track_order):
            row = next(item for item in track_rows if int(item["budget"]) == budget and item["score"] == score)
            value = float(row["precision"])
            y = group_y - track_index * 5
            x_end = e_left + value / 0.75 * e_width
            drawing.add(Line(e_left, y, x_end, y, strokeColor=color, strokeWidth=3))
            text(drawing, x_end + 3, y - 2, f"{value:.2f}", 5.0, DARK)
    for index, (score, color) in enumerate(track_order):
        x = 292 + (index % 2) * 112; y = 25 - (index // 2) * 10
        drawing.add(Line(x, y, x + 10, y, strokeColor=color, strokeWidth=2.5))
        label = {"avi": "AVI", "dnase_only": "DNase", "h3k27ac_only": "H3K27ac", "combined_magnitude": "Combined"}[score]
        text(drawing, x + 14, y - 2.5, label, 5.7, DARK)
    return drawing


def primary_forest_rows(project):
    primary = read_tsv(project / "confirmation/results/primary_by_cell.tsv")
    draws = read_tsv(project / "confirmation/results/bootstrap_draws.tsv.gz")
    one_mb = [row for row in draws if int(row["block_size"]) == 1_000_000]
    rows = []
    for cell in (*CELLS, "confirmation_macro"):
        key = "confirmation_macro_difference" if cell == "confirmation_macro" else f"{cell}_difference"
        values = [float(row[key]) for row in one_mb]
        if cell == "confirmation_macro":
            point = np.mean([float(row["difference"]) for row in primary if row["cell_context"] in CELLS])
        else:
            point = float(next(row["difference"] for row in primary if row["cell_context"] == cell))
        rows.append({
            "cell_context": cell,
            "concordance_minus_magnitude": point,
            "ci95_low": np.quantile(values, 0.025),
            "ci95_high": np.quantile(values, 0.975),
            "bootstrap": "prespecified reranking 1Mb",
        })
    return rows


def figure_3(project, source):
    rank_groups = read_tsv(project / "manuscript/analysis/avi_context_rank_groups.tsv")
    rank_contrasts = read_tsv(project / "manuscript/analysis/avi_context_rank_contrasts.tsv")
    indicators = read_tsv(project / "manuscript/analysis/uncertainty_indicator_models.tsv")
    risk = read_tsv(project / "confirmation/results/risk_coverage.tsv")
    primary = primary_forest_rows(project)

    rank_source = [row for row in rank_groups if int(row["block_size"]) == 1_000_000 and row["cell_context"] in CELLS and row["group"] in ("magnitude_only_top10", "avi_only_top10")]
    rank_contrast_source = [row for row in rank_contrasts if int(row["block_size"]) == 1_000_000]
    agreement_source = []
    for row in indicators:
        if row["indicator"] != "dnase_h3k27ac_sign_disagreement" or row["cell_context"] not in CELLS:
            continue
        agreement_source.append({
            "cell_context": row["cell_context"],
            "agreement_vs_disagreement_odds_ratio": 1 / float(row["odds_ratio"]),
            "ci95_low": 1 / float(row["ci95_high"]),
            "ci95_high": 1 / float(row["ci95_low"]),
            "q_value_bh_20_tests": row["q_value_bh_20_tests"],
        })
    risk_source = []
    for coverage in (0.01, 0.05, 0.10, 0.25, 0.50):
        for policy in ("magnitude", "concordance"):
            values = [float(row[f"{policy}_error"]) for row in risk if row["cell_context"] in CELLS and float(row["coverage"]) == coverage]
            risk_source.append({"coverage": coverage, "policy": policy, "macro_sign_error": np.mean(values)})
    write_tsv(source / "figure_3a_rank_prevalence.tsv", rank_source)
    write_tsv(source / "figure_3b_rank_contrasts.tsv", rank_contrast_source)
    write_tsv(source / "figure_3c_primary_concordance.tsv", primary)
    write_tsv(source / "figure_3d_agreement_marker.tsv", agreement_source)
    write_tsv(source / "figure_3e_selective_risk.tsv", risk_source)

    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Disagreement is informative; concordance leaves ranking unimproved", 13.5, DARK, bold=True)
    panel(drawing, "A", "Top-10% rank discordance separates MPRA yield", 25, 548)
    left, bottom, width, height = 70, 355, 420, 160
    y_max = 0.24
    axes(drawing, left, bottom, width, height, [], [], 0, y_max, [0, 0.05, 0.10, 0.15, 0.20], "emVar prevalence")
    group_width = width / len(CELLS)
    for index, cell in enumerate(CELLS):
        center = left + group_width * (index + 0.5)
        for offset, group, color in ((-10, "magnitude_only_top10", NAVY), (10, "avi_only_top10", ORANGE)):
            row = next(item for item in rank_source if item["cell_context"] == cell and item["group"] == group)
            value = float(row["emvar_prevalence"]); low = float(row["ci95_low"]); high = float(row["ci95_high"])
            x = center + offset; y = bottom + value / y_max * height
            drawing.add(Rect(x - 6, bottom, 12, y - bottom, fillColor=color, strokeColor=None))
            drawing.add(Line(x, bottom + low / y_max * height, x, bottom + high / y_max * height, strokeColor=DARK, strokeWidth=0.8))
        text(drawing, center, bottom - 13, cell, 6.7, DARK, "middle")
    drawing.add(Rect(180, 329, 10, 7, fillColor=NAVY, strokeColor=None)); text(drawing, 195, 329, "Cell-matched only", 6.5)
    drawing.add(Rect(290, 329, 10, 7, fillColor=ORANGE, strokeColor=None)); text(drawing, 305, 329, "AVI only", 6.5)

    panel(drawing, "B", "Cell-matched-only minus AVI-only prevalence", 25, 292)
    rows = [next(row for row in rank_contrast_source if row["cell_context"] == cell) for cell in (*CELLS, "confirmation_macro")]
    forest(drawing, rows, 100, 110, 160, 150, 0.0, 0.18, "prevalence_difference", "difference_ci95_low", "difference_ci95_high", [*CELLS, "Macro"], color=NAVY)
    text(drawing, 180, 88, "Prevalence difference", 7, DARK, "middle")

    panel(drawing, "C", "Prespecified concordance hypothesis", 285, 292)
    forest(drawing, primary, 370, 110, 145, 150, -0.02, 0.025, "concordance_minus_magnitude", "ci95_low", "ci95_high", [*CELLS, "Macro"], color=RED)
    text(drawing, 442, 88, "Sign-error difference", 7, DARK, "middle")

    panel(drawing, "D", "Agreement marks emVar enrichment", 25, 66)
    for index, row in enumerate(agreement_source):
        x = 55 + index * 100
        ratio = float(row["agreement_vs_disagreement_odds_ratio"])
        text(drawing, x, 43, row["cell_context"], 6.5, GRAY, "middle")
        text(drawing, x, 25, f"{ratio:.2f}x", 10, TEAL, "middle", True)
    return drawing


def supplementary_figures(project, source):
    utility = read_tsv(project / "confirmation/results/budget_utility.tsv")
    write_tsv(source / "figure_s1_all_cell_budget.tsv", utility)
    drawings = []
    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Supplementary Figure S1. Frozen policy performance by cell", 13, DARK, bold=True)
    positions = [(25, 325), (285, 325), (25, 70), (285, 70)]
    for panel_index, (cell, (x, y)) in enumerate(zip(CELLS, positions, strict=True)):
        panel(drawing, chr(ord("A") + panel_index), cell, x, y + 220)
        left, bottom, width, height = x + 35, y + 35, 205, 165
        x_values = [left + (math.log10(budget) - math.log10(96)) / (math.log10(5000) - math.log10(96)) * width for budget in BUDGETS]
        axes(drawing, left, bottom, width, height, x_values, [str(b) for b in BUDGETS], 0, 0.95, [0, 0.3, 0.6, 0.9])
        for policy in ("magnitude", "concordance", "avi", "pip", "random"):
            values = [float(next(row["precision"] for row in utility if row["cell_context"] == cell and int(row["budget"]) == budget and row["policy"] == policy)) for budget in BUDGETS]
            points = []
            for px, value in zip(x_values, values, strict=True):
                py = bottom + value / 0.95 * height; points.extend([px, py])
                drawing.add(Circle(px, py, 1.8, fillColor=PALETTE[policy], strokeColor=PALETTE[policy]))
            drawing.add(PolyLine(points, strokeColor=PALETTE[policy], strokeWidth=1.2, fillColor=None))
    drawings.append(("figure_s1", drawing))

    contrasts = read_tsv(project / "manuscript/analysis/budget_precision_contrasts.tsv")
    rank = read_tsv(project / "manuscript/analysis/avi_context_rank_contrasts.tsv")
    primary = json.loads((project / "confirmation/results/confirmation_results.json").read_text())
    sensitivity = []
    for budget in BUDGETS:
        for size in (1_000_000, 5_000_000):
            row = next(item for item in contrasts if item["cell_context"] == "confirmation_macro" and item["comparator"] == "avi" and int(item["budget"]) == budget and int(item["block_size"]) == size)
            sensitivity.append({"analysis": "magnitude_minus_avi", "level": budget, **row})
    for size in (1_000_000, 5_000_000):
        row = next(item for item in rank if item["cell_context"] == "confirmation_macro" and int(item["block_size"]) == size)
        sensitivity.append({"analysis": "magnitude_only_minus_avi_only_top10", "level": "top10", **row})
    sensitivity.extend([
        {"analysis": "primary_concordance_minus_magnitude", "level": "10%", "block_size": 1_000_000, "precision_difference": primary["primary_macro_difference"], "ci95_low": primary["primary_1mb_bootstrap_95_ci"][0], "ci95_high": primary["primary_1mb_bootstrap_95_ci"][1]},
        {"analysis": "primary_concordance_minus_magnitude", "level": "10%", "block_size": 5_000_000, "precision_difference": primary["primary_macro_difference"], "ci95_low": primary["primary_5mb_bootstrap_95_ci"][0], "ci95_high": primary["primary_5mb_bootstrap_95_ci"][1]},
    ])
    write_tsv(source / "figure_s2_block_sensitivity.tsv", sensitivity)
    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Supplementary Figure S2. Dependence sensitivity", 13, DARK, bold=True)
    panel(drawing, "A", "Cell-matched - AVI precision by block size", 25, 545)
    for index, budget in enumerate(BUDGETS):
        y = 490 - index * 80
        text(drawing, 80, y - 2, f"{budget} slots", 7, DARK, "end")
        for offset, size, color in ((7, 1_000_000, NAVY), (-7, 5_000_000, ORANGE)):
            row = next(item for item in contrasts if item["cell_context"] == "confirmation_macro" and item["comparator"] == "avi" and int(item["budget"]) == budget and int(item["block_size"]) == size)
            scale_left, scale_width = 100, 400
            value=float(row["precision_difference"]); low=float(row["ci95_low"]); high=float(row["ci95_high"])
            xv=scale_left+value/0.55*scale_width; xl=scale_left+low/0.55*scale_width; xh=scale_left+high/0.55*scale_width
            drawing.add(Line(xl,y+offset,xh,y+offset,strokeColor=color,strokeWidth=1.2)); drawing.add(Circle(xv,y+offset,2.5,fillColor=color,strokeColor=color))
    drawing.add(Line(100, 165, 500, 165, strokeColor=DARK, strokeWidth=0.7))
    for tick in np.linspace(0, .5, 6):
        x=100+tick/.55*400; drawing.add(Line(x,165,x,162,strokeColor=DARK)); text(drawing,x,151,f"{tick:.1f}",6.5,GRAY,"middle")
    drawing.add(Line(120, 120, 135, 120, strokeColor=NAVY, strokeWidth=2)); text(drawing,140,117,"1 Mb",7)
    drawing.add(Line(220, 120, 235, 120, strokeColor=ORANGE, strokeWidth=2)); text(drawing,240,117,"5 Mb",7)
    text(drawing,300,85,"Nearly identical intervals show robustness to wider genomic dependence.",8,DARK,"middle")
    drawings.append(("figure_s2", drawing))

    indicators = read_tsv(project / "manuscript/analysis/uncertainty_indicator_models.tsv")
    write_tsv(source / "figure_s3_uncertainty_indicators.tsv", indicators)
    drawing = Drawing(WIDTH, HEIGHT)
    text(drawing, 25, 582, "Supplementary Figure S3. Prespecified uncertainty indicators", 11, DARK, bold=True)
    names = [
        "dnase_h3k27ac_sign_disagreement",
        "absolute_quantile_difference_per_sd",
        "dnase_exact_vs_other_direction_disagreement",
        "h3k27ac_exact_vs_other_direction_disagreement",
        "magnitude_only_vs_avi_only_top10",
    ]
    short = ["DNase/H3K27ac sign", "Quantile gap (per SD)", "DNase exact/other sign", "H3K27ac exact/other sign", "Cell-matched-only vs AVI-only"]
    for panel_index, cell in enumerate(CELLS):
        x = 25 if panel_index % 2 == 0 else 285
        y = 320 if panel_index < 2 else 65
        panel(drawing, chr(ord("A") + panel_index), cell, x, y + 220)
        rows = [next(row for row in indicators if row["cell_context"] == cell and row["indicator"] == name) for name in names]
        # Log odds-ratio axis from 0.4 to 5.
        left, bottom, width, height = x + 95, y + 35, 145, 165
        drawing.add(Line(left, bottom, left + width, bottom, strokeColor=DARK, strokeWidth=.7))
        one = left + (math.log(1)-math.log(.4))/(math.log(5)-math.log(.4))*width
        drawing.add(Line(one,bottom,one,bottom+height,strokeColor=GRAY,strokeDashArray=[2,2]))
        gap=height/len(rows)
        for i,(row,label) in enumerate(zip(rows,short,strict=True)):
            py=bottom+height-gap*(i+.5)
            value=float(row["odds_ratio"]); low=float(row["ci95_low"]); high=float(row["ci95_high"])
            conv=lambda v:left+(math.log(v)-math.log(.4))/(math.log(5)-math.log(.4))*width
            drawing.add(Line(conv(low),py,conv(high),py,strokeColor=NAVY,strokeWidth=1)); drawing.add(Circle(conv(value),py,2.2,fillColor=NAVY,strokeColor=NAVY))
            text(drawing,left-5,py-2,label,5.7,DARK,"end")
        for tick in (.5,1,2,4):
            px=left+(math.log(tick)-math.log(.4))/(math.log(5)-math.log(.4))*width
            drawing.add(Line(px,bottom,px,bottom-3,strokeColor=DARK)); text(drawing,px,bottom-12,str(tick),6,GRAY,"middle")
    drawings.append(("figure_s3", drawing))
    return drawings


def main():
    project = Path(__file__).resolve().parents[1]
    manuscript = project / "manuscript"
    figure_dir = manuscript / "figures"
    source_dir = manuscript / "source-data"
    pdf_dir = manuscript / "output/pdf"
    for directory in (figure_dir, source_dir, pdf_dir):
        directory.mkdir(parents=True, exist_ok=True)

    drawings = [
        ("figure_1", figure_1(project)),
        ("figure_2", figure_2(project, source_dir)),
        ("figure_3", figure_3(project, source_dir)),
        *supplementary_figures(project, source_dir),
    ]
    for name, drawing in drawings:
        renderSVG.drawToFile(drawing, str(figure_dir / f"{name}.svg"))

    pdf_path = pdf_dir / "AlphaGenome_Atlas_manuscript_figures.pdf"
    pdf = canvas.Canvas(str(pdf_path), pagesize=(WIDTH, HEIGHT), pageCompression=1)
    for _, drawing in drawings:
        renderPDF.draw(drawing, pdf, 0, 0)
        pdf.showPage()
    pdf.save()

    with tempfile.TemporaryDirectory(prefix="alphagenome-figure-render-") as temporary:
        prefix = Path(temporary) / "page"
        subprocess.run(
            ["pdftoppm", "-png", "-r", "300", str(pdf_path), str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for index, (name, _) in enumerate(drawings, start=1):
            shutil.copyfile(Path(f"{prefix}-{index}.png"), figure_dir / f"{name}.png")

    output_files = sorted(
        [*figure_dir.glob("*.svg"), *figure_dir.glob("*.png"), *source_dir.glob("*.tsv"), pdf_path]
    )
    manifest = {
        "renderer_sha256": sha256(Path(__file__)),
        "analysis_manifest_sha256": sha256(manuscript / "analysis/manifest.json"),
        "figure_count": len(drawings),
        "pdf_pages": len(drawings),
        "files": {str(path.relative_to(manuscript)): sha256(path) for path in output_files},
    }
    (manuscript / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(f"rendered={len(drawings)} pdf={pdf_path}")


if __name__ == "__main__":
    main()
