#!/usr/bin/env python3
"""Build the author-bearing bioRxiv PDF and curated supplementary package."""

from __future__ import annotations

import hashlib
import html
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from pypdf import PdfReader
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
)


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"
OUTPUT = MANUSCRIPT / "submission" / "biorxiv_2026-10-14"
MAIN_STEM = "AlphaGenome_Cell_Context_bioRxiv_PREPRINT_2026-10-14"
SUPP_STEM = "AlphaGenome_Supplementary_Methods_and_Results_2026-10-14"
DATA_STEM = "AlphaGenome_bioRxiv_Supplementary_Data_2026-10-14"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_FILES = {
    "TNR": FONT_DIR / "Times New Roman.ttf",
    "TNR-Bold": FONT_DIR / "Times New Roman Bold.ttf",
    "TNR-Italic": FONT_DIR / "Times New Roman Italic.ttf",
    "TNR-BoldItalic": FONT_DIR / "Times New Roman Bold Italic.ttf",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def register_fonts() -> None:
    for name, path in FONT_FILES.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        pdfmetrics.registerFont(TTFont(name, str(path)))
    pdfmetrics.registerFontFamily(
        "TNR", normal="TNR", bold="TNR-Bold", italic="TNR-Italic", boldItalic="TNR-BoldItalic"
    )


def normalize_dashes(text: str) -> str:
    return (
        text.replace("\u2011", "-")
        .replace("\u2012", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2212", "-")
    )


def inline_markup(text: str) -> str:
    text = normalize_dashes(text)
    code_parts: list[str] = []

    def hold_code(match: re.Match[str]) -> str:
        code_parts.append(html.escape(match.group(1)))
        return f"@@CODE{len(code_parts) - 1}@@"

    text = re.sub(r"`([^`]+)`", hold_code, text)
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    for index, code in enumerate(code_parts):
        text = text.replace(f"@@CODE{index}@@", f'<font name="Courier">{code}</font>')
    return text


def make_styles(supplement: bool = False) -> dict[str, ParagraphStyle]:
    styles = getSampleStyleSheet()
    body_size = 9.7 if supplement else 10.3
    body_leading = 12.8 if supplement else 13.8
    ink = colors.HexColor("#17212B")
    heading = colors.HexColor("#102F3B")
    accent = colors.HexColor("#176F72")
    return {
        "title": ParagraphStyle(
            "PaperTitle",
            parent=styles["Title"],
            fontName="TNR-Bold",
            fontSize=18.2,
            leading=21.2,
            textColor=heading,
            alignment=TA_CENTER,
            spaceAfter=12,
        ),
        "author": ParagraphStyle(
            "Author",
            parent=styles["Normal"],
            fontName="TNR-Bold",
            fontSize=11.2,
            leading=14,
            alignment=TA_CENTER,
            textColor=ink,
            spaceAfter=4,
        ),
        "metadata": ParagraphStyle(
            "Metadata",
            parent=styles["Normal"],
            fontName="TNR",
            fontSize=9.5,
            leading=12.2,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#384650"),
            spaceAfter=4,
        ),
        "h1": ParagraphStyle(
            "Section",
            parent=styles["Heading1"],
            fontName="TNR-Bold",
            fontSize=13.8,
            leading=16.5,
            textColor=heading,
            spaceBefore=14,
            spaceAfter=6,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "Subsection",
            parent=styles["Heading2"],
            fontName="TNR-Bold",
            fontSize=11.2,
            leading=13.5,
            textColor=accent,
            spaceBefore=10,
            spaceAfter=5,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=styles["BodyText"],
            fontName="TNR",
            fontSize=body_size,
            leading=body_leading,
            alignment=TA_LEFT,
            textColor=ink,
            spaceAfter=6.2,
            allowWidows=0,
            allowOrphans=0,
            splitLongWords=True,
        ),
        "formula": ParagraphStyle(
            "Formula",
            parent=styles["BodyText"],
            fontName="Courier",
            fontSize=8.5,
            leading=11,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#184F58"),
            spaceBefore=3,
            spaceAfter=8,
        ),
        "reference": ParagraphStyle(
            "Reference",
            parent=styles["BodyText"],
            fontName="TNR",
            fontSize=8.6,
            leading=10.8,
            leftIndent=16,
            firstLineIndent=-16,
            textColor=ink,
            spaceAfter=4.2,
            splitLongWords=True,
        ),
        "bullet": ParagraphStyle(
            "Bullet",
            parent=styles["BodyText"],
            fontName="TNR",
            fontSize=body_size,
            leading=body_leading,
            leftIndent=16,
            firstLineIndent=-9,
            bulletIndent=3,
            textColor=ink,
            spaceAfter=4,
        ),
        "figure_title": ParagraphStyle(
            "FigureTitle",
            parent=styles["BodyText"],
            fontName="TNR-Bold",
            fontSize=9.2,
            leading=11.2,
            textColor=heading,
            spaceBefore=5,
            spaceAfter=3,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=styles["BodyText"],
            fontName="TNR",
            fontSize=8.25,
            leading=10.1,
            textColor=colors.HexColor("#34444C"),
            spaceAfter=0,
            splitLongWords=True,
        ),
    }


class ManuscriptDocTemplate(BaseDocTemplate):
    def __init__(self, filename: Path, *, title: str, author: str, short_title: str, subject: str):
        super().__init__(
            str(filename),
            pagesize=letter,
            leftMargin=0.72 * inch,
            rightMargin=0.72 * inch,
            topMargin=0.70 * inch,
            bottomMargin=0.68 * inch,
            title=title,
            author=author,
            subject=subject,
        )
        self.short_title = short_title
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="article")
        self.addPageTemplates(PageTemplate(id="article", frames=[frame], onPage=self._decorate_page))

    def _decorate_page(self, canvas, doc) -> None:
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#CFD7DB"))
        canvas.setLineWidth(0.35)
        canvas.line(self.leftMargin, 0.48 * inch, letter[0] - self.rightMargin, 0.48 * inch)
        canvas.setFont("TNR", 7.8)
        canvas.setFillColor(colors.HexColor("#56636B"))
        canvas.drawString(self.leftMargin, 0.32 * inch, self.short_title)
        canvas.drawRightString(letter[0] - self.rightMargin, 0.32 * inch, f"Page {doc.page}")
        canvas.restoreState()


def markdown_story(path: Path, *, supplement: bool = False) -> list:
    styles = make_styles(supplement=supplement)
    story: list = []
    in_references = False
    skip_legends = False
    front_matter = True
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.lower().startswith("## ") and line.lower().endswith("figure legends"):
            skip_legends = True
            continue
        if skip_legends and line.startswith("## References"):
            skip_legends = False
        elif skip_legends:
            continue
        if line.startswith("# "):
            story.append(Paragraph(inline_markup(line[2:]), styles["title"]))
            continue
        if line.startswith("## "):
            front_matter = False
            heading_text = line[3:]
            in_references = heading_text == "References"
            story.append(Paragraph(inline_markup(heading_text), styles["h1"]))
            continue
        if line.startswith("### "):
            story.append(Paragraph(inline_markup(line[4:]), styles["h2"]))
            continue
        if front_matter:
            if line.startswith("**Manuscript status:"):
                continue
            style = styles["author"] if line.startswith("**") else styles["metadata"]
            story.append(Paragraph(inline_markup(line), style))
            continue
        if re.fullmatch(r"`[^`]+`\.?", line):
            story.append(Paragraph(inline_markup(line.strip("`.")), styles["formula"]))
            continue
        if in_references and re.match(r"^\d+\.\s", line):
            story.append(Paragraph(inline_markup(line), styles["reference"]))
            continue
        if line.startswith("- "):
            story.append(Paragraph(inline_markup(line[2:]), styles["bullet"], bulletText="•"))
            continue
        story.append(Paragraph(inline_markup(line), styles["body"]))
    return story


def parse_captions() -> dict[str, tuple[str, str]]:
    captions: dict[str, tuple[str, str]] = {}
    current: str | None = None
    title = ""
    body: list[str] = []
    caption_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (MANUSCRIPT / "MANUSCRIPT.md", MANUSCRIPT / "SUPPLEMENTARY_METHODS_RESULTS.md")
    )
    for raw in caption_text.splitlines():
        match = re.match(r"#{2,3} (Figure \d+|Supplementary Figure S\d+) \| (.+)", raw)
        if match:
            if current is not None:
                captions[current] = (title, " ".join(body).strip())
            current, title, body = match.group(1), match.group(2), []
        elif current is not None and (raw.startswith("# ") or raw.startswith("## ") or raw.startswith("### ")):
            captions[current] = (title, " ".join(body).strip())
            current, title, body = None, "", []
        elif current is not None and raw.strip():
            body.append(raw.strip())
    if current is not None:
        captions[current] = (title, " ".join(body).strip())
    return captions


def add_figure_pages(story: list, labels: list[str], images: list[Path], *, supplement: bool = False) -> None:
    styles = make_styles(supplement=supplement)
    captions = parse_captions()
    max_width = 6.65 * inch
    max_height = 7.25 * inch
    for label, image_path in zip(labels, images, strict=True):
        title, caption = captions[label]
        image = Image(str(image_path))
        scale = min(max_width / image.imageWidth, max_height / image.imageHeight)
        image.drawWidth = image.imageWidth * scale
        image.drawHeight = image.imageHeight * scale
        block = [
            image,
            Spacer(1, 6),
            Paragraph(inline_markup(f"{label}. {title}"), styles["figure_title"]),
            Paragraph(inline_markup(caption), styles["caption"]),
        ]
        story.append(PageBreak())
        story.append(KeepTogether(block))


def build_pdf(
    source: Path,
    output: Path,
    *,
    supplement: bool,
    figure_labels: list[str],
    images: list[Path],
    title: str | None = None,
    short_title: str | None = None,
    subject: str | None = None,
) -> None:
    title = title or (
        "Supplementary methods and results"
        if supplement
        else "Cell-context matching improves AlphaGenome prioritization of experimentally active regulatory variants"
    )
    short_title = short_title or (
        "AlphaGenome supplementary information" if supplement else "Cell-context matching for AlphaGenome prioritization"
    )
    subject = subject or ("Supplementary information" if supplement else "Genomics preprint")
    story = markdown_story(source, supplement=supplement)
    add_figure_pages(story, figure_labels, images, supplement=supplement)
    doc = ManuscriptDocTemplate(
        output,
        title=title,
        author="Olga Lavinda",
        short_title=short_title,
        subject=subject,
    )
    doc.build(story)


TABLE_FILES = {
    "S1_fixed_budget_effect_estimates.tsv": MANUSCRIPT / "analysis/budget_precision_contrasts.tsv",
    "S2_ranking_metrics.tsv": MANUSCRIPT / "analysis/ranking_metrics.tsv",
    "S3a_avi_context_rank_groups.tsv": MANUSCRIPT / "analysis/avi_context_rank_groups.tsv",
    "S3b_avi_context_rank_contrasts.tsv": MANUSCRIPT / "analysis/avi_context_rank_contrasts.tsv",
    "S4a_selective_risk_intervals.tsv": MANUSCRIPT / "analysis/selective_risk_intervals.tsv",
    "S4b_selective_risk_contrasts.tsv": MANUSCRIPT / "analysis/selective_risk_contrasts.tsv",
    "S5_uncertainty_indicator_models.tsv": MANUSCRIPT / "analysis/uncertainty_indicator_models.tsv",
    "S6_budget_utility.tsv": PROJECT / "confirmation/results/budget_utility.tsv",
    "S7a_integrity_audit.json": MANUSCRIPT / "analysis/integrity_audit.json",
    "S7b_analysis_inspection_receipt.json": MANUSCRIPT / "analysis/inspection_receipt.json",
    "S7c_package_inspection_receipt.json": MANUSCRIPT / "PACKAGE_INSPECTION_RECEIPT.json",
    "S8a_wrong_cell_grid.tsv": MANUSCRIPT / "analysis/context-controls/wrong_cell_grid.tsv",
    "S8b_wrong_cell_contrasts.tsv": MANUSCRIPT / "analysis/context-controls/wrong_cell_contrasts.tsv",
    "S9a_single_track_estimates.tsv": MANUSCRIPT / "analysis/context-controls/single_track_estimates.tsv",
    "S9b_single_track_contrasts.tsv": MANUSCRIPT / "analysis/context-controls/single_track_contrasts.tsv",
    "S9c_context_control_verification_receipt.json": MANUSCRIPT / "analysis/context-controls/verification_receipt.json",
}

REPRODUCIBILITY_FILES = {
    "protocol/PROTOCOL_v0.1.md": PROJECT / "PROTOCOL_v0.1.md",
    "protocol/EVIDENCE_LEDGER.md": PROJECT / "EVIDENCE_LEDGER.md",
    "protocol/CONFIRMATION_PREQUERY_FREEZE_2026-09-13.md": PROJECT / "confirmation/CONFIRMATION_PREQUERY_FREEZE_2026-09-13.md",
    "protocol/CONFIRMATION_ANALYSIS_FREEZE_2026-09-13.md": PROJECT / "confirmation/CONFIRMATION_ANALYSIS_FREEZE_2026-09-13.md",
    "results/confirmation_results.json": PROJECT / "confirmation/results/confirmation_results.json",
    "results/result_manifest.json": PROJECT / "confirmation/results/result_manifest.json",
    "results/confirmation_inspection_receipt.json": PROJECT / "confirmation/results/inspection_receipt.json",
    "results/manuscript_analysis_summary.json": MANUSCRIPT / "analysis/manuscript_analysis_summary.json",
    "results/manuscript_analysis_manifest.json": MANUSCRIPT / "analysis/manifest.json",
    "results/CLAIMS_EVIDENCE_TABLE.tsv": MANUSCRIPT / "CLAIMS_EVIDENCE_TABLE.tsv",
    "results/figure_manifest.json": MANUSCRIPT / "figure_manifest.json",
}

CODE_NAMES = (
    "summarize_siraj_mpra.py",
    "liftover_metadata_variants.py",
    "validate_reference_alleles_twobit.py",
    "annotate_noncoding_gencode.py",
    "build_prior_benchmark_exclusions.py",
    "select_confirmation_variants.py",
    "derive_confirmation_prediction_features.py",
    "authorize_confirmation_outcomes.py",
    "extract_confirmation_outcomes.py",
    "analyze_confirmation.py",
    "verify_confirmation_results.py",
    "extract_adastra_exact_cell_evidence.py",
    "analyze_manuscript.py",
    "verify_manuscript_analysis.py",
    "render_manuscript_figures.py",
)


def supplementary_readme() -> str:
    return """# Supplementary Data 1

This archive supports the preprint “Cell-context matching improves AlphaGenome prioritization of experimentally active regulatory variants.”

## Contents

- `Supplementary_Tables/`: machine-readable Supplementary Tables S1-S9. Tables S8-S9 report post hoc controls.

## Access boundary

This is the minimum preprint supplement. The complete reproducibility package is reserved for the journal submission. Raw AlphaGenome Atlas/API outputs, the local prediction cache and credentials are excluded. Source MPRA and ADASTRA data remain available from their cited providers.
"""


def copy_package_tree(root: Path) -> list[Path]:
    written: list[Path] = []
    for relative, source in {f"Supplementary_Tables/{k}": v for k, v in TABLE_FILES.items()}.items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        written.append(destination)
    readme = root / "README.md"
    readme.write_text(supplementary_readme(), encoding="utf-8")
    written.append(readme)
    return written


def scan_public_tree(files: list[Path]) -> None:
    secret = re.compile(r"AIza[0-9A-Za-z_-]{30,}|sk-[0-9A-Za-z_-]{20,}")
    forbidden = ("/Users/", "fussykitty", "/Downloads/", "/Desktop/")
    failures: list[str] = []
    for path in files:
        if path.suffix.lower() not in {".md", ".py", ".json", ".tsv", ".txt"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if secret.search(text) or any(token in text for token in forbidden):
            failures.append(str(path))
    if failures:
        raise RuntimeError("Public-output scan failed: " + ", ".join(failures))


def build_supplementary_archive(output_zip: Path, table_output: Path) -> None:
    table_output.mkdir(parents=True, exist_ok=True)
    for name, source in TABLE_FILES.items():
        shutil.copy2(source, table_output / name)
    with tempfile.TemporaryDirectory(prefix="alphagenome-supplement-") as temp:
        root = Path(temp) / "Supplementary_Data_1"
        files = copy_package_tree(root)
        scan_public_tree(files)
        sums = root / "SHA256SUMS.txt"
        sums.write_text(
            "".join(f"{sha256(path)}  {path.relative_to(root)}\n" for path in sorted(files)),
            encoding="utf-8",
        )
        files.append(sums)
        with zipfile.ZipFile(output_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for path in sorted(files):
                archive.write(path, Path("Supplementary_Data_1") / path.relative_to(root))


def preprint_text(source: str) -> str:
    return source


def main() -> None:
    register_fonts()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    main_pdf = OUTPUT / f"{MAIN_STEM}.pdf"
    supplement_pdf = OUTPUT / f"{SUPP_STEM}.pdf"
    archive = OUTPUT / f"{DATA_STEM}.zip"
    tables = OUTPUT / "Supplementary_Tables_S1-S9"

    main_md = OUTPUT / f"{MAIN_STEM}.md"
    main_md.write_text(preprint_text((MANUSCRIPT / "MANUSCRIPT.md").read_text(encoding="utf-8")), encoding="utf-8")
    build_pdf(
        main_md,
        main_pdf,
        supplement=False,
        figure_labels=["Figure 1", "Figure 2", "Figure 3"],
        images=[MANUSCRIPT / "figures" / f"figure_{index}.png" for index in (1, 2, 3)],
    )
    build_pdf(
        MANUSCRIPT / "SUPPLEMENTARY_METHODS_RESULTS.md",
        supplement_pdf,
        supplement=True,
        figure_labels=["Supplementary Figure S1", "Supplementary Figure S2", "Supplementary Figure S3"],
        images=[MANUSCRIPT / "figures" / f"figure_s{index}.png" for index in (1, 2, 3)],
    )
    build_supplementary_archive(archive, tables)

    shutil.copy2(MANUSCRIPT / "SUPPLEMENTARY_METHODS_RESULTS.md", OUTPUT / f"{SUPP_STEM}.md")
    source_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT, check=True, capture_output=True, text=True
    ).stdout.strip()
    deliverables = [main_pdf, supplement_pdf, archive, *sorted(tables.iterdir())]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit_base": source_commit,
        "prepared_revision_date": "2026-10-14",
        "repository_release_required_before_submission": True,
        "main_pdf_pages": len(PdfReader(str(main_pdf)).pages),
        "supplement_pdf_pages": len(PdfReader(str(supplement_pdf)).pages),
        "author": "Olga Lavinda",
        "biorxiv_category": "Genomics",
        "license_selection": "CC BY 4.0",
        "ai_disclosure_in_main_pdf": False,
        "supplementary_tables": sorted(TABLE_FILES),
        "files": {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in deliverables},
    }
    (OUTPUT / "SUBMISSION_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"main_pdf={main_pdf}")
    print(f"supplement_pdf={supplement_pdf}")
    print(f"supplementary_archive={archive}")
    print(f"supplementary_tables={len(TABLE_FILES)} files")


if __name__ == "__main__":
    main()
