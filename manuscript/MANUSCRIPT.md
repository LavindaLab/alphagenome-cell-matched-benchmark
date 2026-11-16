# Cell Specific AlphaGenome Predictions Outperform the Global Impact Score in a Massively Parallel Reporter Assay Benchmark of Noncoding Variants

**Olga Lavinda**

Department of Chemistry and Biochemistry, Stern College for Women, Yeshiva University, New York, NY 10016, USA

Correspondence: olga.lavinda@yu.edu · ORCID: 0000-0003-3577-5984

## Abstract

Sequence models such as AlphaGenome now predict how noncoding variants alter regulatory activity across the genome, and they offer both a single global impact score (AVI) and separate predictions for individual cell types. A scientist planning a reporter assay in a particular cell line must decide which of these to use when choosing variants to test. We compared the options against a massively parallel reporter assay of 246,674 noncoding single-nucleotide variants in five human cell lines, with the ranking rules, endpoints and thresholds fixed before the held-out outcomes were examined. Ranking variants by the DNase and H3K27ac effects predicted for the assayed cell line gave 64.8% mean precision at 384 assay slots per cell line, compared with 37.6% for AVI (difference, 27.1 percentage points; 95% confidence interval, 22.9–31.2), and returned 417 more expression-modulating variants. In a post hoc control that changed only the cell line of the prediction, predictions for the assayed cell line outperformed predictions for the other cell lines by 17.5 points (15.3–19.7), with a gain that ranged from 4.9 points in SK-N-SH to 37.4 in K562. H3K27ac predictions alone gave similar point estimates to the combination of DNase and H3K27ac. Variants for which the DNase and H3K27ac predictions agreed in sign were 1.4 to 1.8 times more likely to be expression-modulating. Requiring that agreement did not reduce error in predicting the direction of the reporter effect, so the primary hypothesis was not supported. For reporter experiments in these cell lines, cell specific predictions matched to the assay returned 1.7 times as many confirmed regulatory variants per synthesis slot as the global score, and the value of matching depended strongly on the cell line.

## Introduction

Most variants associated with human traits lie in noncoding DNA, where their effects depend on the regulatory state of the cell [3]. Massively parallel reporter assays (MPRAs) test such variants directly by measuring whether each allele changes reporter expression in a chosen cell line. Every assay has a limited number of synthesis slots, so the variants selected for testing determine what the experiment can find.

Sequence models can now guide that selection. AlphaGenome predicts the effect of a variant on many molecular readouts, including chromatin accessibility (DNase) and the active histone mark H3K27ac, separately for individual cell types [1]. The AlphaGenome Atlas provides these predictions precomputed for the whole genome, together with the AlphaGenome Variant Impact score (AVI), a single global score that combines model outputs with conservation and protein-related features [2]. The Atlas presents AVI as the quick route to ranking variants. A scientist designing a reporter assay in one cell line therefore faces a practical choice: rank candidates by the global score, or by the regulatory predictions for the cell line that will be assayed.

The biology suggests the choice should matter. Regulatory elements are active in some cell types and silent in others, so a prediction of accessibility or histone marking in the assayed cell line describes something different from a genome-wide estimate of impact. If cell specific predictions carry information about which variants are active in a given cell line, using the prediction for that cell line should raise the hit rate. If the signal is largely shared across cell types, it should add little. The original AlphaGenome study compared its predictions with MPRA effects [1], and the Atlas study evaluated AVI broadly [2]. We did not find a direct comparison of AVI with cell specific Atlas predictions on the large five cell line MPRA of Siraj and colleagues [3].

We used that MPRA to make the comparison. Throughout, a matched prediction is the AlphaGenome prediction for the cell line in which the reporter assay was performed. Before examining the held-out outcomes we fixed two questions. The first, designated the primary hypothesis, was whether requiring the DNase and H3K27ac predictions to agree in sign would reduce error in predicting the direction of a variant's effect. The second was whether ranking by matched regulatory predictions returns more reporter-active variants than ranking by AVI at a fixed assay budget. Matched predictions returned substantially more hits than AVI at every budget. Sign agreement identified variants more likely to be active, and it did not improve prediction of direction. Two further analyses, added after these results were known and labeled post hoc throughout, asked where the gain came from: one changed only the cell line of the prediction, and the other separated DNase from H3K27ac.

## Results

### Study design

We used the MPRA of Siraj et al., which tested 304,278 variants, comprising fine-mapped variants from complex trait and expression quantitative trait locus studies together with control variants, in K562, HepG2 and SK-N-SH cells, with subsets tested in A549 or HCT116 [3]. We adopted the source authors' definition of an expression-modulating variant (emVar): an element with reporter activity and a significant difference between alleles at a false discovery rate below 0.1. An emVar is a hit in the analyses below.

Genomic blocks of one megabase were divided into a pilot set (10% of blocks) and a held-out set (90%). The ranking rules, cell line assignments, scores, assay budgets and outcome definitions were fixed before the held-out outcomes were examined, and all results reported here come from the held-out set (Fig. 1). After conversion to the current genome build, checking of reference alleles and removal of coding and essential splice site variants, 246,674 noncoding single-nucleotide variants had complete predictions for all five cell lines. The candidates available for ranking were the fine-mapped variants: 176,396 in HepG2, 176,198 in K562, 119,946 in HCT116 and 176,100 in SK-N-SH. These four cell lines met a prespecified minimum number of events and form the confirmatory analysis. A549 did not, and it is reported descriptively.

### Matched cell specific predictions returned more hits than the global score

For each variant and cell line, the matched regulatory score was the larger of the absolute DNase and H3K27ac effect quantiles predicted for that cell line. We ranked each cell line's candidates by this score and, separately, by AVI, and counted the emVars among the top-ranked variants at budgets of 96, 384, 1,000 and 5,000 assay slots per cell line. Neither ranking used MPRA outcomes.

At 384 slots per cell line, the matched ranking reached a mean precision of 64.8% across the four cell lines, compared with 37.6% for AVI, a difference of 27.1 percentage points (95% confidence interval, 22.9–31.2; Table 1, Fig. 2A). In absolute terms it returned 995 emVars where AVI returned 578, a gain of 417 hits spread over 367 additional megabase blocks (Fig. 2D). At 1,000 slots it returned 2,250 compared with 1,277. The matched ranking was ahead at every budget and in every confirmatory cell line, and it also exceeded ranking by fine-mapping posterior probability and random selection (Supplementary Fig. S1; Supplementary Tables S1 and S2). The gap narrowed as the budget grew, from 44.0 points at 96 slots to 12.5 at 5,000, as expected when a ranking exhausts its strongest candidates. Intervals computed with five megabase blocks were nearly identical (Supplementary Fig. S2).

**Table 1 | Mean emVar precision (%) across the four confirmatory cell lines for the matched regulatory score and AVI.** Differences are in percentage points with 95% confidence intervals from a block bootstrap of one megabase blocks.

| Assay slots per cell line | Matched regulatory score | AVI | Difference (95% CI) |
|---|---|---|---|
| 96 | 74.7 | 30.7 | 44.0 (35.5–52.0) |
| 384 | 64.8 | 37.6 | 27.1 (22.9–31.2) |
| 1,000 | 56.3 | 31.9 | 24.3 (21.9–26.8) |
| 5,000 | 35.5 | 22.9 | 12.5 (11.7–13.4) |

For a noncoding reporter experiment in these cell lines, matched DNase and H3K27ac predictions therefore produced about 1.7 times as many confirmed regulatory variants per synthesis slot as AVI. AVI is designed for broader use across coding and noncoding variants, and it remains the natural choice when the relevant cell type is unknown or no suitable prediction is available.

### The gain depended on using the assayed cell line, and its size varied among cell lines

The comparison with AVI changes two things at once: the kind of score and the cell line it refers to. To separate them, we ranked each cell line's candidates five times with the same regulatory score, each time computed from the predictions for a different one of the five cell lines. One of the five rankings is matched and four use predictions for other cell lines. The candidates, budgets, outcome definition and tie-breaking were identical, so only the cell line of the prediction changed. We compared the matched ranking with the mean of the four others. This control was added after the confirmatory analysis was complete and is post hoc (Fig. 2B,C; Supplementary Table S8).

At 384 slots, the matched ranking reached 64.8% mean precision compared with 47.2% for the mean of the other cell lines, a difference of 17.5 points (15.3–19.7). At 1,000 slots the difference was 15.0 points (13.7–16.2; Table 2). Predictions for the other cell lines still outperformed AVI, by 9.6 points at 384 slots and 9.3 at 1,000, so both the kind of score and the choice of cell line contributed to the advantage over AVI.

**Table 2 | Post hoc cell line control.** emVar precision (%) of the ranking that used predictions for the assayed cell line and the mean precision of the four rankings that used predictions for other cell lines. Differences are in percentage points with 95% confidence intervals.

| Assayed cell line | Assay slots | Matched prediction | Other cell lines (mean) | Difference (95% CI) |
|---|---|---|---|---|
| HepG2 | 384 | 63.0 | 54.1 | 8.9 (3.6–14.1) |
| K562 | 384 | 81.0 | 43.6 | 37.4 (32.9–41.9) |
| HCT116 | 384 | 63.0 | 44.1 | 18.9 (15.2–22.8) |
| SK-N-SH | 384 | 52.1 | 47.1 | 4.9 (0.7–9.1) |
| Mean of four cell lines | 384 | 64.8 | 47.2 | 17.5 (15.3–19.7) |
| HepG2 | 1,000 | 56.5 | 46.6 | 9.9 (6.8–13.0) |
| K562 | 1,000 | 72.2 | 40.1 | 32.2 (29.5–34.9) |
| HCT116 | 1,000 | 50.6 | 35.4 | 15.2 (12.9–17.5) |
| SK-N-SH | 1,000 | 45.7 | 42.9 | 2.8 (0.03–5.2) |
| Mean of four cell lines | 1,000 | 56.3 | 41.2 | 15.0 (13.7–16.2) |

The size of the benefit differed widely among cell lines, which is a result in its own right. At 384 slots, matching was worth 37.4 points in K562 and 4.9 points in SK-N-SH, a more than sevenfold range across four well-characterized cell lines scored by the same model with the same rule. The matched ranking was the best of the five in HepG2, K562 and HCT116. In SK-N-SH, predictions for some of the other cell lines did as well as or better than the matched prediction, and at 1,000 slots the lower limit of the interval approached zero. Using the prediction for the assayed cell line helped on average, and how much it helps should be estimated for the cell line of interest.

### H3K27ac predictions carried most of the signal

The regulatory score combines two predictions. A second post hoc control ranked the same candidates by DNase alone and by H3K27ac alone (Fig. 2E; Table 3; Supplementary Table S9). H3K27ac alone reached 66.7% mean precision at 384 slots and 55.9% at 1,000, compared with 64.8% and 56.3% for the combined score. The confidence interval for the difference between the two included zero at both budgets. DNase alone was lower, by 7.6 points at 384 slots and 4.7 at 1,000, and still far above AVI. The pattern was not uniform: H3K27ac alone did particularly well in SK-N-SH, while the combined score was higher in K562 and HCT116 at 1,000 slots. In this benchmark the predicted change in H3K27ac carried most of the information, and combining it with DNase gave no clear further improvement.

**Table 3 | Post hoc readout control.** Mean emVar precision (%) across the four confirmatory cell lines and total emVars returned.

| Ranking | Precision, 384 slots | emVars, 384 slots | Precision, 1,000 slots | emVars, 1,000 slots |
|---|---|---|---|---|
| H3K27ac alone | 66.7 | 1,025 | 55.9 | 2,236 |
| DNase and H3K27ac combined | 64.8 | 995 | 56.3 | 2,250 |
| DNase alone | 57.2 | 878 | 51.6 | 2,062 |
| AVI | 37.6 | 578 | 31.9 | 1,277 |

### Where the two rankings disagreed, the matched ranking was more often right

Variants ranked in the top 10% by one method and outside it by the other show where the methods differ. Averaged across the four cell lines, 14.7% of the variants favored only by the matched ranking were emVars, compared with 6.3% of those favored only by AVI, a difference of 8.4 points (7.8–8.9) and a mean ratio of 2.41 (2.27–2.57; Fig. 3A,B). The difference was present in every confirmatory cell line, with odds ratios from 1.61 to 3.93 (Supplementary Fig. S3; Supplementary Tables S3 and S5). AVI may be responding to conservation, coding or splicing effects that a reporter assay does not measure, so disagreement between the two rankings identifies variants for which they reflect different biology.

### Agreement between DNase and H3K27ac identified active variants and did not improve prediction of direction

The primary hypothesis concerned the direction of a variant's effect: whether the alternative allele raises or lowers reporter expression. Among emVars in elements with enhancer activity, we predicted direction in two ways. The first took the sign of whichever prediction, DNase or H3K27ac, was larger. The second required the two predictions to agree in sign. We compared the error rates of the two approaches among the 10% of variants each ranked most confidently.

Requiring agreement did not reduce error. The mean error rate across the four cell lines was 0.64 percentage points higher with the agreement requirement (95% confidence interval, −0.04 to 1.17; Fig. 3C). The interval includes zero, and the estimate is well below the two point difference that we had specified in advance as practically important. The primary hypothesis was therefore not supported.

Agreement was informative in a different way. Variants for which the DNase and H3K27ac predictions agreed in sign were 1.39 to 1.76 times more likely to be emVars than variants for which they disagreed, with odds ratios from 1.41 to 1.82 (Fig. 3D). In a secondary analysis of allele specific transcription factor binding measured in the same cell lines [7], agreement modestly improved direction accuracy in four cell lines and made little difference in SK-N-SH. That analysis concerns binding in the genome, while the primary endpoint concerns transcription from a reporter. Taken together, agreement between the two predictions marks variants that are more likely to be active, and it is best used as an annotation.

## Discussion

A scientist choosing variants for a reporter assay in a particular cell line will confirm more regulatory variants by ranking with the DNase and H3K27ac predictions for that cell line than by ranking with the global impact score. In this benchmark the difference was large: 24 to 27 percentage points of precision at typical assay budgets, or about 1.7 times as many hits for the same synthesis cost, spread across hundreds of genomic regions. The comparison reflects a decision that users of the Atlas already face, because the global score is the most accessible output and the cell specific predictions require the user to choose a cell type and a readout.

Part of the advantage comes from the cell line itself. When the same score was computed from predictions for each of the five cell lines, the matched prediction outperformed the others by about 17 points on average. The gain was far from uniform, at 37 points in K562 and 5 in SK-N-SH. One possible explanation is that a cell line gains most from a matched prediction when its regulatory landscape is distinct from the other profiled cell lines, and that for SK-N-SH the other cell lines carry much of the same information. These data cannot test that explanation, and identifying when cell type matters most is a clear next question.

The predicted change in H3K27ac carried most of the signal, with DNase adding little on average. This points to the predicted effect on an active histone mark as the most useful single readout for prioritizing variants for an enhancer reporter assay, while the differences among cell lines caution against treating one readout as best everywhere.

The primary hypothesis was not supported, and the result clarifies what agreement between predictions means. Agreement between DNase and H3K27ac identified variants more likely to be active, yet it did not help predict whether the alternative allele would raise or lower reporter expression. Two regulatory predictions pointing the same way is evidence that a variant matters. Carrying that shared direction through a short episomal reporter construct is a further step, and requiring agreement also discards variants with a strong effect on only one readout.

The conclusions are bounded by the assay and the cell lines tested. MPRAs measure activity from an episomal reporter. The candidates were fine-mapped human variants. All five cell lines are extensively profiled and were available to the model during training, so the results describe how to use AlphaGenome well in such cell lines; performance in primary tissues and sparsely profiled cell types remains to be tested. AVI and the cell specific predictions come from the same model, so this study compares outputs within AlphaGenome. The two controls that separate cell line and readout were added after the confirmatory analysis and are post hoc.

AlphaGenome 2 was released on 11 November 2026, after these analyses were completed [9]. Its design conditions predictions on cell identity, and its developers report improved recovery of putatively causal variants in single cell immune data. The benchmark used here does not depend on the model and is released in full, so the same comparison can be applied to the new model. Whether a model conditioned on cell identity recovers cell type specificity where matching was least valuable here, as in SK-N-SH, can now be tested directly.

Predictions should be chosen to match the evidence an experiment will produce. For noncoding reporter assays in these cell lines, regulatory predictions for the assayed cell line substantially increased the number of confirmed variants.

## Methods

### Study design

The protocol was fixed before the held-out outcomes were examined. Genomic blocks of one megabase (hg19 coordinates) were assigned to a pilot set or a held-out set by a deterministic pseudo-random rule based on a SHA-256 hash of the block identifier, which placed 10% of blocks in the pilot set and the remainder in the held-out set. Preliminary runs of 100 and 1,000 variants tested data handling only. A pilot of 10,000 variants drawn from pilot blocks could inform the decision to proceed and could not alter the study population, endpoints, thresholds, ranking rules or analysis. Held-out outcomes were extracted only after the predictions, scoring rules and analysis code were final, and pilot outcomes were not used in the held-out analysis.

### MPRA data and outcomes

MPRA data came from Supplementary Tables 1 and 4 of Siraj et al. [3]; Table 4 contained 304,278 assayed variants. We retained the authors' calls of reporter activity and of emVar status for each cell line, defining an emVar as an active element with significant allelic skew at a false discovery rate below 0.1. Allelic skew was the reporter activity of the alternative allele minus that of the reference allele on a log2 scale. The analysis of direction was limited to emVars in elements with enhancer activity (active elements with a positive activity log2 fold change in that cell line), because an increase in DNase or H3K27ac was not assumed to correspond to increased expression from silencer constructs. Fine-mapped single-nucleotide variants formed the candidates for ranking. Control variants from the source study contributed to estimates of prevalence and were not candidates.

### Coordinates, alleles and definition of noncoding variants

Variants in hg19 coordinates were converted to hg38 with the UCSC chain file hg19ToHg38.over.chain.gz and liftOver [4]. Only variants with a unique mapping and alleles matching the UCSC hg38 reference sequence were retained. We excluded indels, noncanonical alleles, ambiguous or nonprimary mappings, reference mismatches and duplicate representations. Using GENCODE release 50 (GRCh38.p14) [5], a variant was noncoding if it did not overlap the coding sequence of a protein-coding transcript or either intronic nucleotide adjacent to an internal protein-coding exon boundary. We also excluded positions used in earlier AlphaGenome or Atlas benchmarks and published examples. No exclusion used AVI, a regulatory score or an MPRA outcome.

Of 255,148 single-nucleotide variants in held-out blocks, 246,674 were retained. Exclusions were 6,723 coding or essential splice site variants, 1,639 reference mismatches, 83 variants with no unique coordinate conversion, 17 overlaps with earlier benchmarks or examples and 12 nonprimary mappings.

### AlphaGenome Atlas predictions

All predictions were retrieved from the AlphaGenome Atlas on 13 and 14 September 2026 using the AlphaGenome SDK version 0.9.0 (Git commit aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d) [8]. This was before the release of AlphaGenome 2 on 11 November 2026, and all results describe the first AlphaGenome model. Cell lines were matched to Atlas predictions by ontology term: A549, EFO:0001086; HepG2, EFO:0001187; K562, EFO:0002067; HCT116, EFO:0002824; SK-N-SH, EFO:0003072. The corresponding DNase track indices were 107, 109, 118, 146 and 153, and the H3K27ac track indices were 328, 338, 392, 465 and 511. For each variant and cell line we retrieved the DNase and H3K27ac effects from the center mask scorers, their Atlas quantile scores, and AVI. All 246,674 variants returned every required score.

### Scores and ranking rules

For variant v in cell line c, the matched regulatory score was the larger of the two absolute quantile scores:

`M(v,c) = max(|q_DNase(v,c)|, |q_H3K27ac(v,c)|)`

The predicted direction was the sign of the effect for whichever readout gave the larger absolute quantile. The agreement rule required nonzero DNase and H3K27ac effects of the same sign and ranked variants by the smaller of the two absolute quantile scores:

`C(v,c) = min(|q_DNase(v,c)|, |q_H3K27ac(v,c)|)`

AVI ranking used the Atlas AVI quantile. Ties were broken by a deterministic pseudo-random order derived from the variant identifier. Comparison rankings were the fine-mapping posterior inclusion probability (PIP), deterministic pseudo-random selection stratified by chromosome and PIP availability, the agreement rule, and the agreement rule limited to two variants per megabase block. Budgets were 96, 384, 1,000 and 5,000 variants per cell line.

### Endpoints

The primary endpoint was the difference in direction error between the agreement rule and the matched regulatory score, among emVars in elements with enhancer activity, at 10% coverage, averaged with equal weight across the confirmatory cell lines. Negative values favor the agreement rule. Coverage levels of 1%, 5%, 25% and 50% were secondary. A difference of two percentage points was specified in advance as practically important. A549 was designated descriptive by a prespecified rule based on the number of events in the pilot.

The secondary endpoint was emVar precision among the top-ranked candidates at each budget in each cell line. Additional outputs were the number of emVars, the number of megabase blocks containing at least one emVar, average precision and the area under the receiver operating characteristic curve.

Prespecified indicators of uncertainty were sign disagreement between DNase and H3K27ac, the absolute difference between their quantile scores, disagreement in direction between the matched prediction and the median of the predictions for the other cell lines, and discordance between the matched ranking and AVI at the top 10%. No composite score, calibration or threshold was fitted to outcomes.

### Post hoc controls

Two controls were specified on 8 October 2026, after the confirmatory analysis was complete. They left the protocol, endpoints and confirmatory results unchanged.

The cell line control ranked each cell line's candidates five times using the regulatory score computed from the predictions for A549, HepG2, K562, HCT116 or SK-N-SH. The comparator was the unweighted mean precision of the four rankings that used predictions for other cell lines; no single other cell line was selected. The readout control ranked the same candidates by DNase alone and by H3K27ac alone. Both controls used the same candidates, outcome definition, tie-breaking and budgets of 384 and 1,000 variants as the confirmatory analysis, and both reproduced the confirmatory results for the combined score and AVI exactly before the new rankings were evaluated. The released data contain the absolute value of each quantile score, so the readout control addresses precision and hit counts and does not address direction.

### Statistical analysis

We report estimates for each cell line and equal-weight means across the four confirmatory cell lines. Confidence intervals account for the spatial correlation of neighboring variants by resampling genomic blocks within chromosomes (2,000 bootstrap replicates of one megabase blocks, with five megabase blocks as a sensitivity analysis). For the primary endpoint, variants were ranked again within each bootstrap replicate. For precision at fixed budgets and for the comparison of variants favored by only one ranking, intervals were computed for the fixed selected lists.

Indicators of uncertainty were analyzed in separate logistic regressions for each cell line with standard errors clustered by megabase block. The false discovery rate was controlled at 0.05 across the 20 tests in the confirmatory cell lines by the Benjamini–Hochberg method [6]. Interpretation rests on effect sizes and confidence intervals.

### Allele specific transcription factor binding

As a secondary analysis in the genome, we used allele specific transcription factor binding from the ADASTRA resource for the same cell lines, distributed through IGVF [7]. The direction of binding was the sign of the alternative minus the reference effect. A variant was included only when all significant transcription factor observations in that cell line agreed in direction.

### Software

Analyses used Python 3.12.13 and NumPy 2.3.5. All figures derive from the source data tables provided in the repository.

## Data and code availability

The MPRA and ADASTRA data are available from Siraj et al. [3] and IGVF [7]. The protocol, analysis code, derived data, figure source data and Supplementary Tables S1–S9 are available at https://github.com/LavindaLab/alphagenome-cell-matched-benchmark (release v1.1.0-biorxiv). The benchmark dataset is archived at Zenodo (https://doi.org/10.5281/zenodo.23003991). AlphaGenome predictions are available from the provider [8].

## Author information

**Corresponding author:** Olga Lavinda, Department of Chemistry and Biochemistry, Stern College for Women, Yeshiva University, New York, NY 10016, USA; olga.lavinda@yu.edu; ORCID 0000-0003-3577-5984.

## Author contributions

Olga Lavinda: Conceptualization, Data curation, Formal analysis, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing – original draft, and Writing – review and editing.

## Competing interests

The author declares no competing interests.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Acknowledgments

The author thanks the AlphaGenome and AlphaGenome Atlas teams, Siraj and colleagues, and the IGVF and ADASTRA contributors for making the model resource, assay data and validation resources available for research.

## Ethics statement

This computational study used publicly available variant-level assay data and did not analyze participant-level records or identifiable private information. Ethics approval and consent to participate were not applicable to this secondary analysis.

## Figure legends

### Figure 1 | Study design

**A,** Workflow from the MPRA of Siraj et al. through variant filtering and retrieval of AlphaGenome Atlas predictions to the held-out comparison. Ranking rules, endpoints and thresholds were fixed before the held-out outcomes were examined. **B,** Cell lines, the corresponding AlphaGenome ontology term, the role of each cell line in the analysis, the number of candidate variants and the number of emVars with enhancer activity used in the analysis of direction. **C,** Filtering of single-nucleotide variants in the held-out genomic blocks. emVar, expression-modulating variant; AVI, AlphaGenome Variant Impact score.

### Figure 2 | Matched cell specific predictions return more reporter-active variants than the global score

**A,** Mean emVar precision across the four confirmatory cell lines at 96, 384, 1,000 and 5,000 assay slots per cell line for five rankings: the matched regulatory score, the agreement rule (labeled Concordance), AVI, fine-mapping posterior inclusion probability (PIP) and random selection. **B,** Post hoc control. Mean precision of the ranking that used predictions for the assayed cell line and of the mean of the four rankings that used predictions for other cell lines. **C,** Post hoc control. Difference between the matched ranking and the mean of the other cell lines in each cell line, with 95% confidence intervals. **D,** Additional emVars, and additional megabase blocks containing at least one emVar, returned by the matched ranking relative to AVI across the four cell lines. **E,** Post hoc control. Mean precision for rankings by AVI, DNase alone, H3K27ac alone and the combined score.

### Figure 3 | Disagreement between rankings and agreement between predictions

**A,** Proportion of emVars among variants in the top 10% of the matched ranking only and in the top 10% of AVI only. Error bars are 95% confidence intervals. **B,** Difference between the two proportions in each cell line and on average (labeled Macro). **C,** Primary endpoint. Difference in direction error between the agreement rule and the matched regulatory score at 10% coverage. Positive values indicate higher error with the agreement rule. Variants were ranked again within each bootstrap replicate. **D,** Odds ratio for emVar status among variants whose DNase and H3K27ac predictions agree in sign, relative to variants whose predictions disagree.

## References

1. Avsec Z, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature*. 2026;649:1206–1218. doi:10.1038/s41586-025-10014-0.
2. Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind technical report. 2026. Official landing page: https://deepmind.google/science/alphagenome/. Accessed 2026-09-14.
3. Siraj L, Castro RI, Dewey HB, et al. Functional dissection of complex trait variants at single-nucleotide resolution. *Nature*. 2026. doi:10.1038/s41586-026-10121-6.
4. Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. *Genome Research*. 2002;12:996–1006. doi:10.1101/gr.229102.
5. Frankish A, Diekhans M, Ferreira A-M, et al. GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Research*. 2019;47:D766–D773. doi:10.1093/nar/gky955.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289–300.
7. Abramov S, Boytsov A, Bykova D, et al. Landscape of allele-specific transcription factor binding in the human genome. *Nature Communications*. 2021;12:2751. doi:10.1038/s41467-021-23007-0.
8. Google DeepMind. AlphaGenome API README and terms summary. Accessed 2026-09-14. https://github.com/google-deepmind/alphagenome/blob/main/README.md.
9. Google DeepMind. AlphaGenome 2: deciphering how cells read DNA to better understand disease. 11 November 2026. https://deepmind.google/blog/alphagenome-2-deciphering-how-cells-read-dna-to-better-understand-disease/. Accessed 2026-11-16.
