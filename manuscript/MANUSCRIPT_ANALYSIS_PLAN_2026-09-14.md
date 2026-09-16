# Manuscript-quality analysis plan

**Locked:** 2026-09-14, after the frozen confirmation result and substantive-value GO; before running the manuscript analyses below.

## Status and inferential boundary

The confirmatory protocol, primary hypothesis, populations, scores, policies, outcomes, coverage levels, assay budgets and value gate remain immutable. These analyses characterize the already-observed result for reporting; they do not create a second confirmatory test or authorize outcome-aware optimization. The original 1 Mb block-bootstrap primary interval and its 5 Mb sensitivity remain the confirmatory inference.

## Necessary analyses

1. **Fixed-budget effect uncertainty.** For every frozen cell, budget (96, 384, 1,000 and 5,000) and policy, retain the exact selected set. Estimate paired precision contrasts between matched DNase/H3K27ac magnitude and each frozen comparator by chromosome-stratified genomic-block bootstrap, conditional on those fixed selections. Use 2,000 replicates at 1 Mb and repeat at 5 Mb. The central contrast is magnitude minus AVI; all others are supporting baselines.
2. **Practical yield.** Report exact emVar hits, hit precision, cumulative absolute MPRA skew, hit-containing 1 Mb blocks, promoter composition and missed emVars from the frozen utility output. Convert the magnitude-minus-AVI precision contrast into additional discoveries at each fixed budget; do not extrapolate beyond the tested candidate corpus.
3. **Ranking discrimination.** Report AUROC and average precision for magnitude, AVI and PIP on the same eligible fine-mapped cell-specific candidate universe. PIP metrics use observations with non-missing PIP and declare their denominator.
4. **Selective risk.** Attach 1 Mb and 5 Mb block-bootstrap intervals to the frozen magnitude and concordance sign-error curves at 1%, 5%, 10%, 25% and 50% coverage. Preserve the failed 10% concordance hypothesis and the original frozen bootstrap as authoritative.
5. **AVI–context rank disagreement.** Operationalize the prespecified `AVI-high/matched-modality-low` and inverse strata using the already-frozen 10% coverage threshold: both top 10%, magnitude-only top 10%, AVI-only top 10%, or neither. Report emVar prevalence and the paired magnitude-only minus AVI-only contrast with 1 Mb and 5 Mb block-bootstrap intervals in every cell. This is descriptive secondary evidence, not a newly optimized abstention rule.
6. **Named uncertainty indicators.** Evaluate the five prespecified indicators—DNase/H3K27ac sign disagreement, their absolute quantile difference, exact-cell versus median-other-cell DNase direction disagreement, the analogous H3K27ac disagreement, and AVI-versus-magnitude top-10% rank discordance—using cell-specific logistic models with 1 Mb cluster-robust uncertainty. Control Benjamini–Hochberg FDR across the 20 confirmatory-cell tests. A549 remains descriptive.
7. **Integrity audit.** Recheck frozen hashes, gate/confirmation separation, outcome authorization timing, coordinate/reference receipts, variant and variant-cell uniqueness, exact ontology mapping, Atlas completeness, analysis-script identity, and absence of scientific-rule mutation.

## Reporting constraints

- Effect sizes and intervals lead; P values and q values are secondary.
- MPRA is independent episomal assay evidence, not in-vivo truth or clinical validity.
- The fixed selection lists are not retuned inside bootstrap replicates.
- No ER-stress/inflammation, TYR, new model, additional comparator, learned calibrator or alternative outcome is permitted.
- Any additional analysis requires a named manuscript claim that cannot otherwise be evaluated.

