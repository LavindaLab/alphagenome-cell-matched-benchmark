# Protocol v0.1 — external reliability of AlphaGenome Atlas

**Frozen:** 2026-09-13, before any AlphaGenome/Atlas output was joined to MPRA outcomes
**Study type:** retrospective external assay evaluation; noncommercial academic research
**Working title:** *Do AlphaGenome Atlas confidence policies transport? Independent assay reliability and experimental utility across five cellular contexts*
**Status:** frozen for engineering and the 10,000-SNV gate set; amendments require a dated record before confirmation-set access

## Scientific question

When AlphaGenome Atlas assigns a large cell-matched regulatory effect, how reliably does the predicted allele direction and experimental priority transport to an independent MPRA, and does agreement between two mechanistically related Atlas modalities reduce error at the same prediction coverage?

This is not a first abstention analysis. AlphaGenome already reported score thresholds calibrated to eQTL sign accuracy. The new question is external transport across assay and cellular context, followed by decision utility under realistic assay budgets.

## Prespecified hypotheses and estimands

### Primary hypothesis

At 10% coverage, a fixed **concordant DNase–H3K27ac policy** has lower MPRA sign error than a **magnitude-only policy**, macro-averaged across A549, HepG2, K562, HCT116 and SK-N-SH.

The primary estimand is:

`macro mean over cells [sign-error(concordant, 10% coverage) - sign-error(magnitude, 10% coverage)]`

Sign error is evaluated among significant cell-specific emVars in enhancer-active constructs. A negative difference favors concordance. Report the point estimate and a chromosome/block bootstrap 95% confidence interval; do not turn the Atlas quantile into a probability.

### Key secondary estimands

1. Sign error at 1%, 5%, 25% and 50% coverage and the area under the risk–coverage curve from 1% to 50%.
2. emVar precision lift over cell-specific prevalence at each coverage among all measured noncoding SNVs.
3. emVar precision among active elements only, separating allelic-effect detection from element activity.
4. False-negative burden: fraction of observed emVars below each Atlas coverage cutoff.
5. Between-cell heterogeneity of the concordance benefit.
6. Experimental return per slot at budgets of 96, 384, 1,000 and 5,000 within each cell.

### Agent recommendations, not empirical facts

- A clear failure to transport is publishable if confidence intervals exclude an effect large enough to guide experiments.
- “Reliability,” “agreement,” and “assay-conditional evidence” are acceptable terms. “Calibrated probability,” “causal,” “pathogenic,” and “ground truth” are not supported by this design.

## Data and versions

### Primary evidence

Siraj et al., *Nature* (online 2026-02-25), DOI 10.1038/s41586-026-10121-6.

- Supplementary Table 1: fine-mapping/control membership, PIP and credible-set metadata.
- Supplementary Table 4: one row per assayed variant; cell-specific activity, allelic skew, uncertainty, significance and counts.
- Supplementary Table 7: empirical-Bayes activity/skew estimates for sensitivity analysis only.
- GENCODE release 50 (June 2026, GRCh38.p14; Ensembl 116) for the frozen noncoding definition.
- GRCh38 reference sequence matching AlphaGenome coordinates for REF validation.

The Nature article and its supplements are CC BY 4.0. The authors' analysis code repository is GPL-3.0. Exact source URLs, checksums and inspected schemas are recorded in `EVIDENCE_LEDGER.md`.

### Independent endogenous sensitivity evidence

At least one central direction of effect must be checked against a versioned, experimentally supported endogenous tier selected without consulting Atlas results:

1. allele-specific base editing or IGVF Variant-EFFECTS with an allele-level molecular outcome;
2. cell-matched caQTL or allele-specific chromatin/TF binding;
3. element-level CRISPRi for locus support only.

Computational scores, AVI, motif predictions and unvalidated clinical assertions are not independent evidence. The exact IGVF release and counts are frozen at retrieval in an amendment because the catalog is versioned and changes over time.

## Analysis units and eligibility

### Primary unit

One `variant × MPRA cell line` observation. The five exact Atlas ontology terms are:

| MPRA cell | Atlas ontology |
|---|---|
| A549 | EFO:0001086 |
| HepG2 | EFO:0001187 |
| K562 | EFO:0002067 |
| HCT116 | EFO:0002824 |
| SK-N-SH | EFO:0003072 |

### Inclusion

- Canonical biallelic SNV with A/C/G/T REF and ALT.
- GRCh38 coordinate with exact reference-allele match.
- A non-missing Siraj cell-specific activity and skew estimate.
- Noncoding under GENCODE 50: no overlap with a protein-coding CDS and not an essential splice donor/acceptor nucleotide for any protein-coding transcript. UTR, intronic and intergenic variants remain eligible and are reported separately.
- Successful Atlas retrieval for both exact-cell DNase and H3K27ac scorers.

### Exclusion

- Indels or noncanonical alleles; public Atlas coverage is currently SNV-only.
- Any Siraj saturation-mutagenesis follow-up construct rather than the main 304,278-variant screen.
- Variants explicitly used as AlphaGenome or Atlas headline examples/evaluation cases when an exact overlap list can be constructed. Report exclusions by source; do not infer overlap only from publication date.
- REF mismatch, ambiguous liftover, duplicate allele representation or irreconcilable orientation.
- Coding or essential-splice variants under the rule above.
- Perturbed Atlas tracks, unless a separately prespecified sensitivity analysis matches that perturbation.

Do not exclude variants based on AVI, Atlas score, MPRA effect direction or desired result.

## Outcome definitions

The Siraj definitions are retained rather than re-called from raw counts:

- `active`: cell-specific activity call reported by the authors; their method requires activity magnitude and multiple-testing criteria.
- `emVar`: a cell-specific significant allelic-skew call within an active element (FDR threshold 0.1 in the source analysis).
- `log2Skew`: alternative-minus-reference reporter activity. MPRASuite defines it as B minus A, with A labelled reference and B alternative.
- `enhancer-active`: `active = 1` and the reported cell-specific activity log2 fold-change is positive. This restriction is used for the primary directional test because increased DNase/H3K27ac is not expected to share direction with a silencer assay construct.

Primary sign agreement is `sign(Atlas raw effect) = sign(MPRA log2Skew)`. Zero Atlas raw effects are predictions without a direction and count as abstentions, not correct calls. Silencer-active constructs and Table 7 posterior estimates are secondary sensitivity analyses.

## Atlas variables and fixed policies

### Required scorers

For each exact ontology term, retrieve the unperturbed direct-cell:

- DNase CenterMask effect raw score and quantile score;
- H3K27ac histone-ChIP CenterMask effect raw score and quantile score;
- scorer and track metadata sufficient to reproduce the selection.

If more than one qualifying track appears, select deterministically: unperturbed before perturbed, paired-end before single-end when both exist, then lexicographically smallest stable track identifier. Freeze the returned identifiers after the 100-variant metadata check and apply the same identifiers to every batch.

### Magnitude-only policy

For variant `v` in cell `c`:

- magnitude `M = max(abs(q_DNase), abs(q_H3K27ac))`;
- predicted direction is the raw-score sign of the modality attaining `M`;
- ties use DNase.

Select the top `k%` by `M`, breaking ties by SHA-256 of the normalized variant identifier.

### Concordance policy

- eligible to predict only when nonzero DNase and H3K27ac raw scores have the same sign;
- confidence `C = min(abs(q_DNase), abs(q_H3K27ac))`;
- predicted direction is their common sign;
- rank by `C`, with the same deterministic tie-break.

At each requested coverage, take the highest-ranked eligible concordant observations. If concordant observations are fewer than the requested coverage, report attainable coverage and do not backfill with discordant cases.

### Additional uncertainty indicators

Report separately; do not learn or tune a composite:

- DNase–H3K27ac sign disagreement;
- absolute difference between their quantile scores;
- exact-cell versus median other-cell direction disagreement for each modality;
- AVI-high/matched-modality-low and AVI-low/matched-modality-high strata;
- missing secondary modality/track metadata.

No model, learned calibrator or optimized weighting may be trained on Atlas outputs.

## Leakage control and frozen split

1. Construct 1 Mb genomic block identifiers as `chromosome:floor((position-1)/1,000,000)`.
2. Hash each block identifier with the fixed salt `AGATLAS-SIRAJ-2026-V01`.
3. Lowest 10% of hash space forms the **gate/engineering set**; the remaining 90% is the **confirmation set**.
4. The 10,000-SNV pilot is sampled only from gate blocks, stratified by chromosome, cell-measurement availability and Siraj cohort/control category—not by emVar outcome.
5. Gate outcomes may decide whether a full run is worth its query cost, but may not change the primary policies, 10% primary coverage, eligibility rule or primary estimand.
6. Confirmation MPRA outcomes remain unjoined to Atlas scores until a dated run manifest records the scorer IDs, exclusions, code commit and pilot decision.

The paper reports the gate set separately. The final confirmatory inference uses only confirmation blocks.

## Statistical analysis

- Report cell-specific estimates and an equal-weight macro-average across cells; do not let the larger HepG2/K562 cohorts dominate.
- Use 2,000 chromosome-stratified 1 Mb block-bootstrap replicates for confidence intervals. A 5 Mb block bootstrap is the dependency sensitivity analysis.
- For the primary paired policy contrast, bootstrap the difference within identical cell/block samples.
- Test monotonic score–risk association with a prespecified ordinal trend across coverage cutoffs; treat this as secondary.
- Estimate cross-cell heterogeneity with a fixed policy-by-cell interaction and block-robust uncertainty; report effect sizes, not only P values.
- Control the false-discovery rate at 0.05 across the named secondary discordance indicators. Do not correct the single primary comparison.
- Report AUPRC and AUROC for comparability, but do not make either the primary endpoint.
- Report denominators, missingness, attainable coverage and positive prevalence for every estimate.

## Aim 2: experimental-budget utility

Run independent cell-specific campaigns. One slot is one unique variant construct tested in one prespecified cell context. Restrict the candidate universe to eligible fine-mapped SNVs; matched controls estimate null yield but are not selectable candidates.

Evaluate fixed policies at 96, 384, 1,000 and 5,000 slots:

1. deterministic pseudo-random selection, stratified by chromosome and PIP availability;
2. maximum PIP only (maximum across duplicated association records);
3. AVI only;
4. matched DNase/H3K27ac magnitude policy;
5. concordance policy;
6. concordance policy with a cap of two variants per 1 Mb block.

Primary utility outcome is cell-specific emVar hit count. Secondary outcomes are cumulative absolute MPRA skew, number of independent 1 Mb blocks with at least one hit, promoter/distal balance and missed-hit count. Use inverse-probability weights if any policy is evaluated on a stratified pilot rather than the full eligible cohort.

## Query and storage plan

| Stage | Unique Atlas requests | Purpose | Gate |
|---|---:|---|---|
| Metadata check | 100 | coordinate/orientation, scorer IDs, schema, cache receipt | ≥99% valid REF/ALT and complete required scores |
| Engineering | 1,000 | retries, rate, filtered response size, deterministic rerun | ≥99% success after bounded retries; byte-identical normalized rerun |
| Locked pilot | 10,000 | scientific feasibility on gate blocks | rules below |
| Full ceiling | 281,327 | all Siraj SNVs before noncoding/overlap filtering | only after written pilot decision |

Each cache record must include normalized variant, requested ontology terms and scorers, raw and quantile values, track/scorer metadata, request time, package commit/version, API/resource version when exposed, response checksum and retry history. Store no key or token.

The 88.5 GB static AVI archive is not part of the pilot. Download it only if the full-run decision shows that local AVI joining saves more time/API volume than it costs in storage and transfer.

## Pilot continuation rules

Proceed beyond 10,000 requests only if all conditions hold:

1. at least 99% of canonical SNVs reconcile to GRCh38 and required Atlas fields after documented exclusions;
2. at least 99% of requested variants return complete required scores after bounded retries;
3. the gate set contains at least 200 eligible enhancer-active emVars in each evaluated cell, or the affected cell is explicitly downgraded before confirmation access;
4. bootstrap intervals show either a potentially decision-relevant concordance effect (absolute sign-error change at least 2 percentage points remains plausible) or a sufficiently precise null/failure result to support a transportability paper;
5. a full run is expected to narrow a named scientific uncertainty, rather than merely add rows;
6. no new literature result occupies the external five-cell transportability and budget-utility claim.

Stop or reduce scope if the API cannot provide stable scorer identifiers, exact context matching fails, the MPRA orientation cannot be reproduced, or the full run would violate rate/term constraints. A null biological result is not itself a stop condition.

## Licensing, conflicts and publication

- All raw Atlas work is noncommercial academic research under the current terms.
- Health AI receives only published aggregate findings/methods or permissive-use AVI-derived work. Academic raw Atlas data, caches and feature importances are not used in commercial deliverables.
- Disclose the investigator's Health AI ownership and any related decision-audit/experiment-prioritization service.
- Do not make clinical, diagnostic or ACMG claims.
- Refresh literature, terms, software versions and publication status immediately before confirmation access and immediately before preprint posting.

## Amendments

An amendment may add a sensitivity analysis or document a source version. It may not silently change the primary population, outcome, policies, coverage, estimand or split after gate outcomes are inspected. Every amendment must state what was known, why the change was necessary and whether confirmation data had been accessed.
