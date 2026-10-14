# Cell-context matching AlphaGenome benchmark

Public reproducibility package for:

> Olga Lavinda. *Cell-context matching improves AlphaGenome prioritization of experimentally active regulatory variants.*

This repository preserves the frozen protocol, outcome-access records, analysis code, row-level derived analysis tables, aggregate results, figure source data, and claim-to-artifact mapping used by the manuscript.

## Primary entry points

- `SCREENING_README.md`: concise screener-facing overview and reproduction guide.
- `manuscript/MANUSCRIPT.md`: author manuscript source.
- `manuscript/submission/biorxiv_2026-10-14/AlphaGenome_Cell_Context_bioRxiv_PREPRINT_2026-10-14.pdf`: revised submission PDF.
- `manuscript/CLAIMS_EVIDENCE_TABLE.tsv`: claim-level evidence map.
- `protocol/PROTOCOL_v0.1.md`: frozen study protocol.
- `confirmation/CONFIRMATION_PREQUERY_FREEZE_2026-09-13.md`: confirmation query freeze.
- `confirmation/CONFIRMATION_ANALYSIS_FREEZE_2026-09-13.md`: confirmation analysis freeze.
- `confirmation/prediction_features.tsv.gz`: row-level derived AlphaGenome ranking features.
- `confirmation/avi_scores.tsv.gz`: row-level AVI scores.
- `confirmation/mpra_outcomes_confirmation.tsv.gz`: confirmation-set MPRA outcomes derived from the cited source dataset.
- `confirmation/results/confirmation_results.json`: frozen confirmatory result record.
- `manuscript/analysis/`: manuscript-level analyses and integrity audit.
- `manuscript/analysis/context-controls/`: post hoc exact-cell, nonmatching-context and single-track controls.
- `tools/`: retrieval, processing, analysis, verification, and figure scripts.

## Reproduction boundary

The repository includes the row-level derived data needed to inspect and rerun the reported analysis. Raw AlphaGenome API response objects and the local response cache are excluded. Repeating retrieval from the provider requires independent AlphaGenome access and compliance with the current AlphaGenome Terms of Service.

The source MPRA and ADASTRA records remain attributable to their cited providers. See `DATA_USE_NOTICE.md` for the precise source and reuse boundary.

## Verification

The manuscript analysis can be rerun from the repository root with:

```bash
python3 tools/analyze_manuscript.py --project-dir . --output-dir manuscript/analysis-reproduced
python3 tools/verify_confirmation_results.py \
  confirmation/results \
  confirmation/OUTCOME_ACCESS_AUTHORIZATION.json \
  confirmation/results/inspection_receipt.reproduced.json
```

The post hoc context controls can be checked with:

```bash
python3 tools/verify_context_controls.py
```

The earlier submitted version is tagged `v1.0.1-biorxiv`. The revised package is versioned as `v1.1.0-biorxiv` and adds the exact-cell versus nonmatching-context control, the DNase-only and H3K27ac-only controls, and an updated Figure 2.

## External source records

- Siraj MPRA article and supplementary data: https://doi.org/10.1038/s41586-026-10121-6
- IGVF ADASTRA release: https://data.igvf.org/reference-files/IGVFFI5943XCOS/
- AlphaGenome software and access documentation: https://github.com/google-deepmind/alphagenome
