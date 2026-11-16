# Overview of the study

This study asks a practical question for anyone designing a reporter assay in a particular cell line: should candidate variants be ranked by the AlphaGenome global impact score (AVI), or by the AlphaGenome regulatory predictions for the cell line that will be assayed?

## What was evaluated

The outcome is the massively parallel reporter assay of Siraj et al. An expression-modulating variant (emVar) is an element with reporter activity and a significant difference between alleles at the source study's false discovery rate threshold. The analysis covers 246,674 noncoding single-nucleotide variants with predictions for five cell lines. HepG2, K562, HCT116 and SK-N-SH form the confirmatory analysis, and A549 is reported descriptively.

Two questions were fixed before the held-out outcomes were examined.

The first, designated primary, was whether requiring the DNase and H3K27ac predictions to agree in sign reduces error in predicting the direction of a variant's effect. It did not: at 10% coverage the mean error rate was 0.64 percentage points higher with the agreement requirement (95% confidence interval, -0.04 to 1.17). Agreement did identify variants more likely to be active.

The second was whether ranking by the predictions for the assayed cell line returns more reporter-active variants than ranking by AVI at a fixed assay budget. It did, in every confirmatory cell line:

| Assay slots per cell line | Matched predictions | AVI | Difference (95% CI) | Additional emVars |
| ---: | ---: | ---: | ---: | ---: |
| 384 | 64.8% | 37.6% | 27.1 points (22.9 to 31.2) | 417 |
| 1,000 | 56.3% | 31.9% | 24.3 points (21.9 to 26.8) | 973 |

Two controls were added after these results were known and are labeled post hoc. The first ranked the same candidates with the same score computed from the predictions for each of the five cell lines. Predictions for the assayed cell line outperformed the mean of the other cell lines by 17.5 points at 384 slots (15.3 to 19.7) and 15.0 points at 1,000 (13.7 to 16.2), with the smallest gain in SK-N-SH. The second ranked by DNase alone and by H3K27ac alone: H3K27ac alone gave similar point estimates to the combined score, and DNase alone was lower.

## What can be rerun from this repository

The repository contains the frozen protocol, analysis code, row-level derived
prediction and outcome tables, aggregate result tables, and verification
receipts needed to rerun the reported analysis locally. The analysis command
does not retrieve AlphaGenome data or require the excluded raw API response
cache. Repeating provider retrieval is outside this package and requires
independent AlphaGenome access under the provider's current terms. Source MPRA
and ADASTRA records remain subject to their cited providers' reuse terms; see
`DATA_USE_NOTICE.md`.

## Run from the repository root

```bash
python3 tools/analyze_manuscript.py --project-dir . --output-dir manuscript/analysis-reproduced
python3 tools/verify_confirmation_results.py \
  confirmation/results \
  confirmation/OUTCOME_ACCESS_AUTHORIZATION.json \
  confirmation/results/inspection_receipt.reproduced.json
```

The post hoc controls are verified with:

```bash
python3 tools/verify_context_controls.py
```

Expected analysis output begins with:

```text
passed budget384=0.271484 ci=[0.229340,0.312221] budget1000=0.243250 ci=[0.218678,0.268218]
```

The confirmation verifier should report `passed=True checks=25/25`, and the
context-control verifier should report `passed=True checks=15/15`. On the
validation machine used for this package, the manuscript analysis completed in
approximately 96 seconds with the committed derived tables; runtime will vary
with hardware and filesystem speed.
