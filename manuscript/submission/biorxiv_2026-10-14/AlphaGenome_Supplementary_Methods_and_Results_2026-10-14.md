# Supplementary methods and results

## Supplementary methods

### Frozen roles and populations

The gate-pilot continuation rule required at least 200 eligible enhancer-active emVars per evaluated cell or a pre-confirmation downgrade. A549 did not meet this rule in the gate analysis and was therefore fixed as descriptive before confirmation outcomes were opened. HepG2, K562, HCT116 and SK-N-SH form every confirmatory macro estimate with equal weight.

The measured, noncoding confirmation populations contained 245,897 HepG2, 245,606 K562, 151,396 HCT116 and 245,542 SK-N-SH variant–cell observations. Their fine-mapped candidate universes contained 176,396, 176,198, 119,946 and 176,100 variants, respectively. The directional enhancer-active emVar populations were 13,045, 15,809, 5,866 and 13,762. A549 counts were 100,611 measured noncoding observations, 62,319 fine-mapped candidates and 2,249 directional observations.

### Deterministic policy selection

Every policy operated independently within a cell and selected exactly the requested budget where available. PIP was reduced to the maximum value across duplicated association records. Random selection was deterministic and stratified by chromosome and PIP availability. Concordance did not backfill with discordant variants. The block-capped concordance sensitivity retained at most two selected variants in each frozen hg19 1 Mb block. All ties used the same hash of the normalized hg38 variant identifier.

### Bootstrap interpretation

The primary directional bootstrap rebuilt cell-specific rankings after sampling genomic blocks and is the only confirmatory inference for the concordance hypothesis. The later manuscript bootstrap held selected lists fixed, because its purpose was to attach sampling uncertainty to the already-frozen experimental allocations. The two procedures answer different questions. A conditional 10%-coverage interval from the later fixed-list risk analysis excludes zero, but it must not replace the original reranking interval, which includes zero. The manuscript therefore reports the concordance hypothesis as failed.

### Cluster-robust indicator models

Each named indicator was entered alone into a logistic regression with an intercept and emVar status as the response. The sandwich covariance clustered observations by frozen hg19 1 Mb block. Absolute DNase–H3K27ac quantile difference was standardized within cell. The 20 tests comprised five indicators in each of four confirmatory cells; A549 q values were not used for confirmatory interpretation. All 20 confirmatory associations remained below FDR 0.05; effect sizes and consistency governed interpretation.

### Post hoc context and single-track controls

The context and single-track controls were specified on 8 October 2026 after completion of the confirmatory record and initial manuscript. They reused the released prediction features and held-out MPRA outcomes without changing eligibility, outcome definitions, budgets or tie-breaking. The context control ranked each outcome cell's fixed candidate set using combined DNase/H3K27ac magnitude predictions from all five available cell contexts. Its fixed comparator was the unweighted mean of the four nonmatching-context precisions. The single-track control reconstructed absolute DNase and H3K27ac quantiles from the released maximum, minimum and winning-modality fields, then ranked the same candidates by DNase alone, H3K27ac alone, combined magnitude and AVI. Reconstruction identities passed for all 1,233,370 prediction rows at an absolute tolerance of 1e-12.

Both controls used 2,000 chromosome-stratified block-bootstrap replicates conditional on each fixed ranked selection at 1 Mb and 5 Mb. The 1 Mb and 5 Mb seeds were 21,261,008 and 25,261,008, respectively. These analyses are post hoc and address prioritization precision and hit counts. Signed single-track quantiles were unavailable in the released feature representation, so single-track directional accuracy was not evaluated.

## Supplementary results

### Budget robustness and frozen baselines

The equal-weight magnitude-minus-AVI precision differences were 44.01 percentage points at 96 slots (1 Mb 95% CI, 35.50–52.02), 27.15 at 384 (22.93–31.22), 24.32 at 1,000 (21.87–26.82) and 12.53 at 5,000 (11.70–13.36). The 5 Mb intervals were 35.37–51.80, 22.99–31.44, 21.96–26.94 and 11.66–13.43, respectively.

Magnitude ranking also exceeded PIP and random selection at every budget. Against concordance, magnitude had 5.21 points higher precision at 96 slots (0.13–10.40), 0.91 at 384 (−1.55 to 3.27), 1.52 at 1,000 (0.14–2.84) and 2.45 at 5,000 (1.95–2.91). These hit-precision comparisons do not change the negative directional concordance result.

### Genome-wide ranking discrimination

Matched magnitude exceeded AVI in both average precision and AUROC in all confirmatory cells. Average precision was 0.214 versus 0.139 in HepG2, 0.279 versus 0.137 in K562, 0.207 versus 0.113 in HCT116 and 0.171 versus 0.132 in SK-N-SH. AUROC was 0.714 versus 0.651, 0.729 versus 0.632, 0.730 versus 0.662 and 0.677 versus 0.633. PIP average precision ranged from 0.068 to 0.088 and AUROC from 0.553 to 0.598 across these cells.

### Selective risk

At the prespecified 10% coverage, the frozen reranking primary contrast was +0.64 percentage points (1 Mb 95% CI, −0.04 to 1.17; 5 Mb, −0.004 to 1.12), where positive values favor magnitude. Fixed-list characterization across 1%, 5%, 10%, 25% and 50% coverage produced small, nonmonotonic contrasts, including a negative point difference at 25%. This pattern does not support a portable concordance-based abstention policy.

### Rank-disagreement sensitivity

The equal-weight mean of the four cell-specific cell-matched-only versus AVI-only emVar prevalence ratios was 2.41 (1 Mb 95% CI, 2.27–2.57), with ratios of 2.10 in HepG2, 3.32 in K562, 2.70 in HCT116 and 1.54 in SK-N-SH. Using 5 Mb blocks produced the same point estimates and materially unchanged intervals. This analysis used the already-prespecified 10% threshold and did not search alternative cutoffs.

### Endogenous, mechanistically adjacent sensitivity

The ADASTRA exact-cell consensus set contained 283, 662, 1,031, 213 and 93 variant–cell observations in A549, HepG2, K562, HCT116 and SK-N-SH. Magnitude-policy direction accuracy was 0.625, 0.637, 0.646, 0.587 and 0.710. Among concordance-eligible observations, accuracy was 0.677, 0.657, 0.676, 0.632 and 0.708. These small-sample results are compatible with endpoint-specific effects but do not validate MPRA direction.

### Cellular context and single-track attribution

Exact-cell combined-magnitude ranking exceeded the fixed nonmatching-context mean at both tested budgets. Four-cell macro precision was 0.6478 compared with 0.4723 at 384 slots, a difference of 0.1755 (1 Mb 95% CI, 0.1535–0.1967; 5 Mb, 0.1536–0.1979). At 1,000 slots, precision was 0.5625 compared with 0.4124, a difference of 0.1501 (1 Mb, 0.1375–0.1621; 5 Mb, 0.1373–0.1632). Differences were positive in HepG2, K562, HCT116 and SK-N-SH. Exact-cell ranking was highest among the five prediction contexts in HepG2, K562 and HCT116; in SK-N-SH, some individual nonmatching contexts equaled or exceeded it while their fixed mean remained lower.

H3K27ac-only macro precision was 0.6673 at 384 slots and 0.5590 at 1,000, compared with 0.6478 and 0.5625 for the combined score. H3K27ac-minus-combined differences were +0.0195 (1 Mb 95% CI, −0.0041 to 0.0442) and −0.0035 (−0.0157 to 0.0082), respectively. DNase-only macro precision was 0.5716 and 0.5155, lower than the combined score by 0.0762 (−0.0963 to −0.0567) and 0.0470 (−0.0586 to −0.0358). The cell-level pattern was heterogeneous, supporting H3K27ac as the dominant macro signal in this benchmark without establishing one universally optimal track.

## Supplementary tables

- Supplementary Table S1: all fixed-budget effect estimates and 1 Mb/5 Mb intervals (`analysis/budget_precision_contrasts.tsv`).
- Supplementary Table S2: ranking metrics (`analysis/ranking_metrics.tsv`).
- Supplementary Table S3: top-10% disagreement groups and contrasts (`analysis/avi_context_rank_groups.tsv`, `analysis/avi_context_rank_contrasts.tsv`).
- Supplementary Table S4: selective-risk intervals and contrasts (`analysis/selective_risk_intervals.tsv`, `analysis/selective_risk_contrasts.tsv`).
- Supplementary Table S5: cluster-robust uncertainty-indicator models (`analysis/uncertainty_indicator_models.tsv`).
- Supplementary Table S6: frozen exact-budget utility outcomes (`confirmation/results/budget_utility.tsv`).
- Supplementary Table S7: integrity audit and inspection receipts (`analysis/integrity_audit.json`, `analysis/inspection_receipt.json`).
- Supplementary Table S8: post hoc exact-cell versus nonmatching-context estimates and contrasts (`analysis/context-controls/wrong_cell_grid.tsv`, `analysis/context-controls/wrong_cell_contrasts.tsv`).
- Supplementary Table S9: post hoc single-track estimates and contrasts (`analysis/context-controls/single_track_estimates.tsv`, `analysis/context-controls/single_track_contrasts.tsv`).

## Supplementary figure legends

### Supplementary Figure S1 | Fixed-policy precision by outcome cell

Cell-specific emVar precision for exact-cell magnitude, concordance, AVI, PIP and deterministic pseudo-random ranking at 96, 384, 1,000 and 5,000 assay slots. The same fixed policies and candidate definitions were applied within each cell.

### Supplementary Figure S2 | Sensitivity to genomic dependence block size

Equal-weight four-cell magnitude-minus-AVI precision differences with chromosome-stratified block-bootstrap intervals using 1 Mb and 5 Mb blocks. Similar intervals across block sizes indicate that the fixed-budget comparison is robust to the wider dependence assumption.

### Supplementary Figure S3 | Prespecified uncertainty indicators

Cell-specific logistic-regression odds ratios and 95% confidence intervals for emVar status. Standard errors were clustered by fixed hg19 1 Mb block, and Benjamini–Hochberg correction was applied across the 20 confirmatory-cell tests. The top-10% comparison contrasts variants selected only by cell-matched magnitude with those selected only by AVI.
