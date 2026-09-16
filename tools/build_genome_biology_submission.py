#!/usr/bin/env python3
"""Build the Genome Biology review manuscript from the verified core source."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from build_biorxiv_submission import MANUSCRIPT, build_pdf, register_fonts, sha256


OUTPUT = MANUSCRIPT / "submission" / "genome_biology_2026-09-14"
STEM = "AlphaGenome_Cell_Matched_Genome_Biology_SUBMISSION_2026-09-14"
TITLE = "Cell-matched AlphaGenome predictions improve experimental variant prioritization over a global impact score"


def journal_text(source: str) -> str:
    declaration = """## Declaration of generative AI and AI-assisted technologies in the writing process

During preparation of this manuscript, OpenAI Codex assisted with literature organization, code review, computational verification, figure production and language editing. The author reviewed and revised all outputs, verified reported results against the frozen records and takes full responsibility for the content.

"""
    marker = "## Figure legends\n"
    if marker not in source:
        raise RuntimeError("Figure-legends marker is missing")
    return source.replace(marker, declaration + marker, 1)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    register_fonts()
    source_path = OUTPUT / f"{STEM}.md"
    pdf_path = OUTPUT / f"{STEM}.pdf"
    source_path.write_text(journal_text((MANUSCRIPT / "MANUSCRIPT.md").read_text(encoding="utf-8")), encoding="utf-8")
    build_pdf(
        source_path,
        pdf_path,
        supplement=False,
        figure_labels=["Figure 1", "Figure 2", "Figure 3"],
        images=[MANUSCRIPT / "figures" / f"figure_{index}.png" for index in (1, 2, 3)],
        title=TITLE,
        short_title="Cell-matched AlphaGenome prioritization",
        subject="Genome Biology submission manuscript",
    )
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "journal": "Genome Biology",
        "article_type": "Research",
        "scientific_source": str(MANUSCRIPT / "MANUSCRIPT.md"),
        "journal_disclosure_added": True,
        "files": {
            source_path.name: {"bytes": source_path.stat().st_size, "sha256": sha256(source_path)},
            pdf_path.name: {"bytes": pdf_path.stat().st_size, "sha256": sha256(pdf_path)},
        },
    }
    (OUTPUT / "SUBMISSION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(pdf_path)


if __name__ == "__main__":
    main()
