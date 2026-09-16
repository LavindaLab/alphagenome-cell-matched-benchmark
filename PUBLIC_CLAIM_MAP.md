# Public claim-to-artifact map

The frozen, row-level claim map is `manuscript/CLAIMS_EVIDENCE_TABLE.tsv`. Canonical paths in that table are preserved in this repository.

The principal records are:

| Manuscript result | Evidence record | Analysis or verification code |
|---|---|---|
| Negative prespecified concordance hypothesis | `confirmation/results/confirmation_results.json`; `confirmation/results/bootstrap_draws.tsv.gz` | `tools/analyze_confirmation.py`; `tools/verify_confirmation_results.py` |
| Cell-matched versus AVI fixed-budget results | `manuscript/analysis/budget_precision_contrasts.tsv`; `confirmation/results/budget_utility.tsv` | `tools/analyze_manuscript.py`; `tools/verify_manuscript_analysis.py` |
| Rank-disagreement analysis | `manuscript/analysis/avi_context_rank_groups.tsv`; `manuscript/analysis/avi_context_rank_contrasts.tsv` | `tools/analyze_manuscript.py` |
| Row-level prediction and outcome inputs | `confirmation/prediction_features.tsv.gz`; `confirmation/avi_scores.tsv.gz`; `confirmation/mpra_outcomes_confirmation.tsv.gz` | retrieval and processing scripts under `tools/` |
| Retrieval completeness and integrity | `confirmation/atlas/retrieval_summary.json`; `manuscript/analysis/integrity_audit.json` | `tools/verify_confirmation_retrieval.py`; `tools/verify_manuscript_analysis.py` |

