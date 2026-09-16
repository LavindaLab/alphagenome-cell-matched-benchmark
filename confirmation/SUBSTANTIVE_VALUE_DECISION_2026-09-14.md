# Substantive-value decision

**Decision:** **GO — narrowly scoped academic manuscript expansion is scientifically justified.**  
**Decision date:** 2026-09-14  
**Scope completed:** frozen held-out confirmation and the pre-specified post-result value gate  
**Manuscript status:** not started by this gate

This is a GO because candidate finding A passes its frozen quantitative screen and all five mandatory qualitative requirements. Candidate finding C independently passes its quantitative screen and strengthens the usage rule. The frozen primary concordance hypothesis (candidate B) does not pass and must be reported as such.

## What the held-out evaluation established

### A. Practical experimental-prioritization result — passes

Among noncoding variants eligible for each cell, selecting the largest exact-cell DNase/H3K27ac AlphaGenome effect materially increased emVar return over global AVI at both frozen assay budgets.

| Budget per cell | Magnitude macro precision | AVI macro precision | Absolute gain | Relative gain | Confirmatory cells favoring magnitude |
|---:|---:|---:|---:|---:|---:|
| 384 | 64.78% | 37.63% | 27.15 pp | 72.1% | 4/4 |
| 1,000 | 56.25% | 31.93% | 24.32 pp | 76.2% | 4/4 |

At 384 slots per cell, magnitude selection produced 995 emVar hits versus 578 for AVI: **417 additional hits**. Those hits occupied 896 versus 529 independent 1 Mb blocks: **367 additional hit-containing blocks**. At 1,000 slots, the totals were 2,250 versus 1,277 hits (**973 additional**) across 1,764 versus 1,090 hit-containing blocks (**674 additional**). The effect is therefore not a duplicate-selection artifact or a one-cell result.

The concrete usage change is: for a cell-matched noncoding reporter assay with a limited synthesis budget and relevant exact-cell tracks, rank first by the largest matched DNase/H3K27ac effect rather than by global AVI alone. AVI remains appropriate as a broad genome-wide triage score when the assay context is unknown, an exact track is unavailable, or coding/general deleteriousness is the target.

### B. Frozen primary concordance hypothesis — fails

The primary four-cell macro difference in sign error was `sign-error(concordance) - sign-error(magnitude) = +0.006374`, so concordance was 0.64 percentage points worse, not better. Its 1 Mb block-bootstrap 95% interval was `[-0.000399, 0.011682]`; the 5 Mb sensitivity interval was `[-0.000037, 0.011153]`. Both include zero and the absolute point estimate is below the frozen 2 percentage-point importance threshold. A claim that cross-modal concordance improves high-confidence MPRA direction calls is unsupported.

### C. Agreement-marker/ranking-policy distinction — passes

DNase/H3K27ac sign agreement was associated with higher emVar prevalence by 1.70-fold in HCT116, 1.73-fold in HepG2, 1.76-fold in K562 and 1.39-fold in SK-N-SH. Thus the frozen 1.5-fold threshold was met in 3/4 confirmatory cells. At the same time, the primary 1 Mb interval rules out a concordance benefit of 2 percentage points or more; the 5 Mb sensitivity reaches the same practical conclusion.

Agreement is therefore useful as a coarse reliability annotation, but replacing magnitude ranking with a concordance-only high-confidence policy did not improve reporter-assay direction. Exact-cell ADASTRA allele-specific TF-binding evidence is compatible with this endpoint-specific interpretation: concordance modestly improved binding-direction accuracy in four cell lines and was essentially unchanged in SK-N-SH, without establishing that TF-binding concordance must improve episomal transcriptional direction.

## Mandatory qualitative gate

| Requirement | Decision | Evidence and boundary |
|---|---|---|
| Effect magnitude | Pass | Candidate A exceeds the pre-specified relative and absolute thresholds at both budgets in all 4/4 confirmatory cells, with hundreds of additional hits and independent loci. Candidate C meets its prevalence threshold in 3/4 large cell-specific populations. |
| Novelty | Pass, medium risk | The 2026-09-14 refresh distinguishes fixed-budget Siraj/exact-cell/AVI utility from AlphaGenome CAGI5 MPRA, Atlas AVI/SGE, RNU4-2, Weykopf and Nadal Vidal. Broad “first benchmark,” “first MPRA,” “first multimodal” and universal-calibration claims are prohibited. |
| Experimental utility | Pass | It changes which score should fill a 384- or 1,000-construct cell-matched assay and preserves substantially more locus diversity. |
| Interpretability | Pass | MPRA regulatory activity is mechanistically closer to exact-cell accessibility/active-histone effects than to an 18-feature, cross-context deleteriousness score. The no-material-concordance-benefit conclusion survives 5 Mb resampling and is not contradicted by the adjacent ADASTRA sensitivity. This is an assay-conditional explanation, not proof of an in-vivo causal mechanism. |
| Change to use | Pass | Prefer matched molecular-effect magnitude for context-specific noncoding assay allocation; treat agreement as annotation, not the ranking rule; retain AVI for broad or context-unknown triage; do not automate when the relevant context/track is absent. |

## What changed from the original brief

- Generic abstention/calibration was dropped because AlphaGenome already reports confidence/coverage behavior and adjacent calibration work exists.
- The proposed concordance policy is not the positive headline; the held-out result rejects its expected practical advantage.
- The publishable contribution is the large fixed-budget advantage of exact-cell molecular effects over global AVI, plus the marker-versus-ranking distinction.
- MPRA remains independent episomal evidence, not in-vivo ground truth. ADASTRA remains a mechanistically adjacent sensitivity, not a substitute.

## Provenance and integrity

- 246,674 frozen variants were retrieved with 246,674 complete records, zero incomplete records and two retries.
- The run made 246,674 Atlas base requests at a 900/minute cap, took 16,450 seconds (4 h 34 min) for retrieval, and produced a 91,282,147-byte compressed raw cache (87.1 MiB).
- Prediction derivation produced 1,233,370 variant-cell rows and did not access MPRA outcomes.
- Outcome access occurred only after the prediction freeze and checksum authorization.
- The frozen analysis used confirmation outcomes but permanently excluded gate outcomes; scientific rules were not mutated.
- Independent result verification passed 25/25 checks. Result commit: `0422dbf7678f0c84fda9d6a9336fc7b1b10ad8e2`.
- Frozen protocol SHA-256: `09bcfa0c09d6c39f6b98c521d73e772629395252f6faa4a35e4d3a8314d43f5f`.
- Frozen pre-outcome evidence-ledger SHA-256: `05df7123ab873d25b3d4b93d0e93dcfdeccad98b009c8fe55e73b39dff7f382c`.
- Frozen analysis-script SHA-256: `f34b02fea83d59360ea4efadccfe9018c37a5910b8c033f1f69820bc2f03440e`.

## Academic and Health AI lanes

**Academic:** GO for a concise external decision-utility paper or research note. Lead with assay-slot return and the context-matched-versus-global-score decision; report the failed primary hypothesis prominently.

**Health AI:** HOLD raw Atlas/API data, cached outputs and derived non-permissive features out of commercial systems. Current official documentation limits API and Atlas information to noncommercial use except named permissive downloadable artifacts and prohibits using outputs to train another ML model. A future commercial product may use the published aggregate decision rule, separately licensed model access or permissive inputs only after a fresh terms/licensing review. This is an operational boundary, not legal advice.

## Next scientific gate

Before public posting, independently reproduce the manuscript-facing tables/figures from the committed result files, verify every reported number and caption against the frozen outputs, refresh forward citations once more, and confirm the release package excludes raw/noncommercial Atlas outputs and secrets. No new outcome-aware optimization is justified.

