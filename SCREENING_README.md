# AlphaGenome cell-context screening package

This package tests a biological design question: when the experiment is a
cell-specific MPRA, does matching AlphaGenome regulatory predictions to the
assayed cell line improve the recovery of expression-modulating variants? A
prespecified comparison evaluates exact-cell DNase/H3K27ac ranking against the
global AlphaGenome Variant Impact score (AVI). Post hoc controls hold the score
fixed while changing prediction context and separate DNase from H3K27ac.

## What was evaluated

The independent outcome is the Siraj et al. MPRA: an expression-modulating
variant (emVar) is an active element with significant allelic skew at the
source study's FDR threshold. The confirmation analysis covers 246,674
noncoding SNVs with complete derived predictions across five matched cell
contexts. HepG2, K562, HCT116, and SK-N-SH are the prespecified confirmatory
cells; A549 is descriptive under the prespecified event-count rule.

The prespecified primary endpoint failed. At 10% coverage among enhancer-active emVars,
sign-concordance had a four-cell macro sign-error difference of +0.00637
relative to magnitude ranking (concordance minus magnitude; 1-Mb block-bootstrap
95% CI, −0.00040 to 0.01168). The interval includes zero and the point
estimate is below the prespecified practical-importance threshold. This package
therefore does not support a general concordance advantage for MPRA direction.

The prespecified secondary endpoint is experimental yield under fixed
screening budgets. Cell-matched magnitude ranking outperformed AVI in every
confirmatory context:

| Budget per cell | Cell-matched precision | AVI precision | Difference (95% CI) | Additional emVars |
| ---: | ---: | ---: | ---: | ---: |
| 384 | 0.647786 | 0.376302 | 0.271484 (0.229340, 0.312221) | 417 |
| 1,000 | 0.562500 | 0.319250 | 0.243250 (0.218678, 0.268218) | 973 |

The post hoc context control ranked the same candidate sets with the same
combined score computed in five cellular contexts. Exact-cell mean precision
was 0.6478 versus 0.4723 for the fixed mean of four nonmatching contexts at 384
slots (difference 0.1755; 1 Mb 95% CI, 0.1535–0.1967) and 0.5625 versus 0.4124
at 1,000 slots (difference 0.1501; 0.1375–0.1621). The difference was positive
in all four outcome cells, with the smallest effect in SK-N-SH. A post hoc
single-track control found similar macro performance for H3K27ac alone and the
combined score at both budgets, while DNase alone was lower. These controls do
not alter the prespecified primary result or secondary AVI comparison.

## Exact reproduction boundary

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
