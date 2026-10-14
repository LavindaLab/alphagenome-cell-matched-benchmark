# Cell-context matching improves AlphaGenome prioritization of experimentally active regulatory variants

**Olga Lavinda**

Department of Chemistry and Biochemistry, Stern College for Women, Yeshiva University, New York, NY 10016, USA

Correspondence: olga.lavinda@yu.edu · ORCID: 0000-0003-3577-5984

## Abstract

AlphaGenome Atlas offers several predictions per variant, creating a choice for cell-specific experiments. We evaluated that choice against an independent massively parallel reporter assay (MPRA) covering 246,674 noncoding single-nucleotide variants in five human cell lines, with policies, endpoints and thresholds fixed before confirmatory outcomes were accessed. In a post hoc control that held the candidate set and scoring rule constant while changing only prediction context, exact-cell ranking reached 64.78% mean precision at 384 assay slots compared with 47.23% for the unweighted mean of four nonmatching contexts, a difference of 17.55 percentage points (1 Mb block-bootstrap 95% confidence interval, 15.35–19.67). At 1,000 slots the difference was 15.01 points (13.75–16.21). The benefit was positive in all four confirmatory cells, although it was small in SK-N-SH and some individual nonmatching contexts equaled or exceeded the exact context there. In the prespecified secondary comparison, exact-cell ranking reached 64.78% precision compared with 37.63% for the global AlphaGenome Variant Impact score (AVI), returning 417 additional expression-modulating variants across 367 additional megabase blocks at 384 slots and 973 across 674 blocks at 1,000. A post hoc single-track analysis found similar macro precision for H3K27ac alone and the combined DNase/H3K27ac score, while DNase alone was lower. The prespecified primary hypothesis failed: requiring DNase and H3K27ac effects to agree in sign made MPRA direction error 0.64 percentage points worse than magnitude ranking (−0.04 to 1.17). Across the tested cell lines, matching regulatory predictions to the assayed context improved experimental yield; cross-modal agreement annotated regulatory plausibility without improving directional ranking.

## Introduction

AlphaGenome Atlas supplies molecular effects for particular cell types and assays alongside AVI, a global score combining AlphaGenome outputs with conservation and protein-related features [1,2]. For a scientist allocating a cell-specific reporter assay, these outputs create a biological design question: which prediction should determine the variants that are tested? Genome-wide precomputation removes much of the inference burden from variant screening, while the Atlas interface foregrounds AVI as the single-number route to rapid ranking [2]. Whether a broad score or a biologically matched output better guides a context-specific experiment remains an empirical question.

Regulatory elements are context-specific. A prediction of accessibility or active-histone marking in the assayed cell line therefore represents a different biological proposition from a genome-wide impact estimate. If context-specific predictions carry information about which variants will be active in a given cellular environment, matching them to the assay should raise the hit rate. If the signal is largely shared across contexts, matching should add little.

AVI supports broad prioritization across coding and noncoding mechanisms. Exact-cell DNase and H3K27ac predictions represent regulatory processes closer to the evidence produced by an MPRA. AlphaGenome's original study benchmarked MPRA effects, including matched DNase predictions and learned combinations [1]; the Atlas study benchmarked AVI broadly and against saturation genome-editing screens [2]. The documented literature search through 2026-09-14 did not identify a zero-shot comparison of AVI and exact-cell Atlas outputs on the large five-cell MPRA reported by Siraj and colleagues [3].

We used that MPRA for an outcome-blind evaluation. The prespecified primary hypothesis asked whether requiring DNase and H3K27ac effects to agree in sign would reduce directional error at fixed coverage; a prespecified secondary aim compared experimental yield under fixed ranking policies and assay budgets. The primary hypothesis failed. Matched ranking nevertheless produced substantially more MPRA hits than AVI. Post hoc controls then tested whether this gain depended on cellular context and whether it required combining DNase with H3K27ac. Together, the analyses support matching the prediction to the assayed cell, choosing a regulatory output aligned with the experiment, and using cross-modal agreement as an annotation.

## Results

### An independent reporter assay provided held-out outcomes

We used the 304,278-variant MPRA of Siraj et al., which tested fine-mapped trait and expression quantitative trait locus variants and matched controls in K562, HepG2 and SK-N-SH, with subsets in A549 or HCT116 [3]. We retained the source authors’ cell-specific definition of an expression-modulating variant (emVar): an active element with significant allelic skew at false-discovery rate below 0.1.

Policy definitions, contexts, scores, coverage thresholds, assay budgets and outcome definitions were fixed before confirmation-outcome access (Fig. 1). Genomic blocks were assigned deterministically to a 10% gate set or 90% confirmation set, and gate outcomes informed only the decision to run the held-out analysis. A prespecified event-count rule retained HepG2, K562, HCT116 and SK-N-SH for equal-weight confirmatory summaries; A549 remained descriptive.

After coordinate conversion, reference-allele checking and removal of coding and essential-splice variants, 246,674 unique noncoding SNVs carried complete predictions across all five exact-cell ontologies (1,233,370 variant–cell rows). Independent inspections passed every retrieval, confirmation-result and manuscript-analysis check; no duplicate row or gate/confirmation crossover was detected. Eligible fine-mapped candidate universes contained 176,396 variants for HepG2, 176,198 for K562, 119,946 for HCT116 and 176,100 for SK-N-SH.

### Matching the prediction to the assayed cell line raised the hit rate

To isolate cellular context from other properties of the score, we ranked each cell's fixed candidate set five times using the same combined DNase/H3K27ac magnitude score computed in five prediction contexts. We compared the exact-cell ranking with the unweighted mean of the four nonmatching contexts. Candidate sets, eligibility rules, outcome definitions, budgets and tie-breaking were identical across the five rankings; only the cellular context of the prediction changed. No nonmatching context was selected on the basis of outcomes. This control was specified after the confirmatory record closed and is post hoc throughout (Fig. 2B,C; Supplementary Table S8).

Exact-cell ranking reached a four-cell macro precision of 0.6478 at 384 assay slots compared with 0.4723 for the fixed all-wrong-cell mean, a difference of 0.1755 (1 Mb block-bootstrap 95% CI, 0.1535–0.1967; 5 Mb, 0.1536–0.1979). At 1,000 slots, precision was 0.5625 compared with 0.4124, a difference of 0.1501 (0.1375–0.1621; 5 Mb, 0.1373–0.1632).

The difference was positive in every confirmatory cell and varied substantially in size. At 384 slots the cell-level differences were 0.0892 for HepG2 (0.0355–0.1409), 0.3737 for K562 (0.3293–0.4190), 0.1895 for HCT116 (0.1525–0.2276) and 0.0495 for SK-N-SH (0.0067–0.0908); at 1,000 slots they were 0.0988 (0.0681–0.1303), 0.3215 (0.2949–0.3490), 0.1518 (0.1289–0.1748) and 0.0283 (0.0003–0.0521). Exact-cell ranking was highest among the five prediction contexts in HepG2, K562 and HCT116. In SK-N-SH, some individual nonmatching contexts matched or exceeded the exact context even though their fixed mean was lower. The control therefore supports an average benefit of exact matching across these four outcome cells and five tested contexts.

### Matched regulatory ranking outperformed the global impact score

The cell-matched magnitude policy ranked variants by the larger absolute quantile score from the exact-cell DNase and H3K27ac predictions. AVI ranked the same eligible fine-mapped, noncoding candidates by the Atlas global-score quantile. Neither policy used MPRA outcomes for construction or threshold selection.

At 384 selected variants per cell, cell-matched magnitude achieved a four-cell macro precision of 0.6478 versus 0.3763 for AVI, an absolute difference of 0.2715 (1 Mb block-bootstrap 95% CI, 0.2293–0.3122; Fig. 2A). At 1,000 variants, precision was 0.5625 versus 0.3193, a difference of 0.2432 (0.2187–0.2682; Fig. 2A). The advantage was positive in every confirmatory context (Supplementary Table S1).

Cell-matched ranking remained better than AVI, the fixed PIP baseline and random selection at all four budgets (96, 384, 1,000 and 5,000; Supplementary Fig. S1 and Supplementary Table S1). The macro precision difference over AVI declined from 0.4401 at 96 slots (0.3550–0.5202) to 0.1253 at 5,000 (0.1170–0.1336). AVI precision rose between 96 and 384 slots because ranks 97–384 had higher aggregate emVar prevalence than its top 96. Five-megabase dependence blocks gave nearly identical intervals (Supplementary Fig. S2). Average precision and AUROC also favored cell-matched ranking in each cell (Supplementary Table S2).

The AVI comparison changes both prediction context and score composition. The post hoc control separates the context component while keeping the combined regulatory score fixed. At 384 slots, the exact-versus-wrong-context contrast accounted arithmetically for 0.1755 of the 0.2715 point-estimate difference between exact-cell ranking and AVI; the remaining 0.0960 separated wrong-context regulatory ranking from AVI. At 1,000 slots the corresponding point-estimate components were 0.1501 and 0.0932. Confidence intervals were computed for the exact-versus-wrong-context component, not for this two-part arithmetic decomposition.

At 384 slots per cell, cell-matched ranking returned 995 emVars compared with 578 for AVI: 417 additional hits across 367 additional hit-containing 1 Mb blocks (Fig. 2D). At 1,000 slots it returned 2,250 compared with 1,277, a gain of 973 across 674 additional hit-containing blocks. The additional experimental yield therefore extended across hundreds of genomic blocks.

For a noncoding reporter experiment in these contexts, exact-cell accessibility and active-histone effects produced more confirmed regulatory variants per synthesis slot than AVI. AVI addresses broader coding and noncoding impact and remains relevant when the context is unknown, a suitable track is unavailable or the intended evidence is outside regulatory activity.

### H3K27ac carried most of the matched-score signal

The ranking policy combined two modalities by taking the larger absolute quantile. A post hoc single-track control, specified with the context control, ranked the same candidates by DNase alone and H3K27ac alone (Fig. 2E; Supplementary Table S9). Absolute single-track quantiles were reconstructed from the released fields for all 1,233,370 prediction rows with zero identity failures at a tolerance of 1 × 10⁻¹².

H3K27ac alone reached a macro precision of 0.6673 at 384 slots and 0.5590 at 1,000, compared with 0.6478 and 0.5625 for the combined score. The H3K27ac-minus-combined difference was +0.0195 at 384 (1 Mb, −0.0041 to 0.0442; 5 Mb, −0.0034 to 0.0413), corresponding to 30 more hits, and −0.0035 at 1,000 (1 Mb, −0.0157 to 0.0082; 5 Mb, −0.0146 to 0.0087), corresponding to 14 fewer. DNase alone retained prioritization signal at a lower level, 0.0762 below the combined score at 384 (1 Mb, −0.0963 to −0.0567) and 0.0470 below at 1,000 (1 Mb, −0.0586 to −0.0358). The cell-level pattern was heterogeneous: H3K27ac-only was particularly favorable in SK-N-SH, while the combined score was higher in K562 and HCT116 at 1,000 slots.

The combined maximum showed no clear macro improvement over H3K27ac alone at either budget. Thus the macro gain reflects exact-cell regulatory prediction, with H3K27ac capturing nearly all of the observed prioritization signal; the data do not identify a general advantage from taking the maximum across the two modalities.

### Rank disagreement locates where the global score misses

At the prespecified 10% threshold, we classified variants as highly ranked by both policies, cell-matched only, AVI only or neither. This descriptive comparison used fixed scores and did not optimize a threshold. Across confirmatory contexts, equal-weight mean emVar prevalence was 0.1471 among cell-matched-only variants and 0.0633 among AVI-only variants, a difference of 0.0838 (1 Mb block-bootstrap 95% CI, 0.0784–0.0891). The equal-weight mean of the four cell-specific prevalence ratios was 2.41 (2.27–2.57; Fig. 3A,B); it is distinct from the ratio of the two macro prevalences.

The separation was present in every confirmatory context; cluster-robust odds ratios for cell-matched-only versus AVI-only membership ranged from 1.61 to 3.93 after Benjamini–Hochberg control across the 20 prespecified indicator tests (Supplementary Fig. S3 and Supplementary Tables S3, S5).

DNase–H3K27ac sign disagreement, a larger cross-modal quantile gap and exact-cell versus other-cell direction disagreement were also associated with lower emVar odds in all four contexts. In this dataset, disagreement functions as a warning annotation; calibration and portable thresholds remain open questions.

### The prespecified concordance hypothesis failed

The primary endpoint was MPRA allelic-direction error at 10% coverage among enhancer-active emVars. Magnitude predicted direction from the modality with the larger effect. Concordance required nonzero DNase and H3K27ac effects with the same sign and ranked eligible observations by the smaller absolute quantile score.

The hypothesis failed. Four-cell macro sign error was 0.64 percentage points worse under concordance than magnitude ranking (prespecified contrast, +0.00637; reranking 1 Mb block-bootstrap 95% CI, −0.00040 to 0.01168; Fig. 3C). The 5 Mb interval was −0.00004 to 0.01115. Both intervals included zero, and the point estimate was well below the prespecified 2-point threshold for practical importance. Cell-specific estimates had the same direction, while fixed-list risk–coverage differences were small, coverage-dependent and changed sign (Supplementary Table S4). The MPRA results provide no evidence of a general concordance advantage.

Agreement marked a different property: emVar prevalence was 1.39–1.76 times higher among sign-agreeing than sign-disagreeing variants, with corresponding logistic-regression odds ratios of 1.41–1.82 (Fig. 3D). In exact-cell ADASTRA allele-specific transcription-factor binding, concordance modestly improved direction accuracy in four lines and was nearly neutral in SK-N-SH. That evidence concerns endogenous binding, while the primary endpoint concerns reporter transcription. Agreement therefore annotates regulatory plausibility without improving MPRA directional ranking.

## Discussion

The prespecified primary hypothesis failed. Requiring DNase and H3K27ac effects to agree in sign enriched for emVars and left reporter allelic-direction error unimproved relative to magnitude ranking. Agreement records whether two regulatory layers point in the same direction; carrying that shared sign through a 200-bp episomal reporter is a further step the data did not support. Requiring agreement also excludes strong single-modality effects. The small ADASTRA sensitivity analysis was compatible with modest directional information for endogenous binding, and its different endpoint leaves the MPRA result intact. In this study, agreement is an annotation.

The supported biological finding concerns cellular context. When a fixed candidate set was ranked with the same regulatory score computed in five prediction contexts, exact-cell ranking exceeded the mean of the four nonmatching contexts by 17.6 precision points at 384 slots and 15.0 at 1,000. The difference was positive in every confirmatory cell. Its size varied markedly, and individual nonmatching contexts sometimes matched the exact context in SK-N-SH. The evidence therefore supports an average benefit across these four outcome cells and five prediction contexts, with weaker context specificity in SK-N-SH.

Against the global impact score, matched ranking improved hit precision by 24–27 percentage points at two assay budgets and returned hundreds of additional hits across hundreds of additional genomic blocks. This comparison addresses the practical choice between an accessible global score and outputs aligned with a regulatory assay. Because prediction modality and cellular context change together in that comparison, the wrong-cell control is essential: it shows that context matching contributes substantially to the observed difference without treating the point-estimate decomposition as causal.

The single-track control sharpens the biological interpretation. H3K27ac captured nearly all of the combined score's macro prioritization performance, while DNase carried a weaker signal. Taking the maximum across both tracks gave no clear macro improvement over H3K27ac alone. This result points to matched active-histone prediction as the dominant signal in this MPRA benchmark, while the heterogeneous cell-level results caution against treating one track as universally optimal.

Rank disagreement makes the assay alignment operational. Variants selected only by matched regulatory magnitude had more than twice the emVar prevalence of variants selected only by AVI at the same 10% threshold. AVI may capture conservation, coding, splicing or cross-context effects outside the MPRA endpoint, so disagreement identifies candidates for which the two evidence layers prioritize different biology.

The inference is bounded by the assay, candidate set and tested contexts. MPRA measurements are episomal and assay-conditional, supplying independent experimental evidence about reporter activity. Candidates were enriched for fine-mapped human variants and controls. All five cell lines are extensively profiled lines with available DNase and H3K27ac tracks; transfer to primary tissues, perturbations and sparsely profiled contexts remains to be tested. AVI and the molecular scores come from the same model resource, so this study compares decision layers within AlphaGenome. Fixed-budget intervals condition on fixed selected lists; the primary directional interval reranked within bootstrap replicates and remains authoritative. The context and single-track controls were specified after the confirmatory record closed and retain post hoc status. The documented literature search ran through 2026-09-14.

Prediction layers should be chosen for the evidence an experiment will produce. In these cell-specific noncoding reporter assays, matching regulatory predictions to the assayed context changed experimental yield substantially. A plausible cross-modal concordance requirement left directional reliability unimproved.

## Methods

### Study design and outcome firewall

The protocol was frozen before confirmation-outcome access. One-megabase hg19 blocks, defined as `chromosome:floor((position−1)/1,000,000)`, were hashed with SHA-256 and salt `AGATLAS-SIRAJ-2026-V01`. Digests in the lowest 10% of hash space formed the gate set; all others formed confirmation. The 100- and 1,000-variant stages tested only schema, coordinates, alleles, context, retrieval, caching and computational plumbing. A 10,000-SNV gate-block pilot could inform the run/no-run decision but not the primary population, endpoint, thresholds, policies or confirmation analysis. Confirmation outcomes were extracted only after scorer identities, predictions, code hashes and analysis rules were frozen and inspected; gate outcomes were never imported.

### MPRA evidence and outcomes

External evidence came from Supplementary Tables 1 and 4 of Siraj et al. [3]; Table 4 contained 304,278 assayed variants. We retained the authors’ cell-specific `active` and emVar calls, defining an emVar as an active element with significant allelic skew at FDR <0.1. Following MPRASuite, `Log2Skew` was alternative-minus-reference reporter activity. The directional endpoint was limited to enhancer-active emVars (`active=1` and positive cell-specific activity log2 fold change), because increased DNase/H3K27ac direction was not assumed to map to silencer constructs. Eligible fine-mapped SNVs formed the prioritization population; matched controls informed prevalence but were not selectable.

### Coordinates, alleles and noncoding definition

Source hg19 SNVs were represented as half-open one-base intervals and converted to hg38 with the UCSC `hg19ToHg38.over.chain.gz` chain and command-line liftOver [4]. Only unique mappings with alleles matching the official UCSC hg38 twoBit sequence were retained. We excluded noncanonical alleles, indels, nonprimary or ambiguous mappings, reference mismatches and duplicate allele representations. Under GENCODE release 50 (GRCh38.p14/Ensembl 116), noncoding variants did not overlap a protein-coding transcript CDS or either intronic nucleotide adjacent to an internal protein-coding exon boundary [5]. We also excluded exact coordinates used in located AlphaGenome/Atlas benchmarks or headline examples. No exclusion used AVI, a molecular score or an MPRA outcome.

Of 255,148 canonical SNVs in confirmation blocks, 246,674 remained. The frozen pipeline excluded 83 unmapped or nonunique lifts, 1,639 hg38 reference mismatches, 6,723 CDS/essential-splice variants, 12 nonprimary mappings and 17 prior-benchmark/example overlaps.

### AlphaGenome Atlas retrieval and context matching

We used the official AlphaGenome SDK version 0.9.0 at Git commit `aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d`. Exact ontology terms were A549 EFO:0001086, HepG2 EFO:0001187, K562 EFO:0002067, HCT116 EFO:0002824 and SK-N-SH EFO:0003072. Frozen DNase track indices were 107, 109, 118, 146 and 153 in that order; H3K27ac indices were 328, 338, 392, 465 and 511. Required scorers were unperturbed direct-cell DNase CenterMask effect and H3K27ac histone-ChIP CenterMask effect, with deterministic tie-breaking among qualifying tracks. Each request also returned AVI. We retained raw effects, Atlas quantiles and track metadata.

Retrieval used four workers, at most 900 requests per minute and two bounded retries. It made 246,674 base requests, required two retries, completed in 16,450 seconds and produced a 91,282,147-byte compressed research cache. All required scores were finite. API credentials came from the environment and were not persisted.

### Frozen scores and policies

For variant `v` in cell `c`, matched magnitude was

`M(v,c) = max(|q_DNase(v,c)|, |q_H3K27ac(v,c)|)`.

Direction was the sign of the maximizing modality’s raw effect, with DNase breaking exact ties. Concordance required nonzero DNase and H3K27ac effects with the same sign and used

`C(v,c) = min(|q_DNase(v,c)|, |q_H3K27ac(v,c)|)`.

AVI ranking used the Atlas AVI quantile. SHA-256 of the normalized variant identifier broke ties. Frozen baselines were chromosome- and PIP-availability-stratified deterministic pseudo-random selection, maximum PIP, concordance, and concordance capped at two variants per 1 Mb block. Cell-specific budgets were 96, 384, 1,000 and 5,000 variants.

### Endpoints

The primary estimand was the equal-weight mean across confirmatory cells of `sign-error(concordance) − sign-error(magnitude)` at 10% coverage among enhancer-active emVars; negative values favored concordance. A549 was descriptive under the frozen pilot count rule. Secondary coverage levels were 1%, 5%, 25% and 50%.

Budget utility was evaluated on each cell’s eligible fine-mapped population. The principal measure was emVar precision; secondary outputs were emVar count, cumulative absolute MPRA skew, number of 1 Mb hit-containing blocks, promoter composition and missed hits. Promoter-proximal meant within 2 kb of a GENCODE 50 protein-coding transcript start site.

Prespecified indicators were DNase–H3K27ac sign disagreement, absolute cross-modal quantile difference, exact-cell versus median-other-cell direction disagreement for each modality, and AVI/matched top-10% rank discordance. We fit no composite, calibrator or outcome-trained threshold.

### Post hoc context and single-track controls

Two diagnostics were specified on 8 October 2026, after the confirmatory endpoint, result record and initial manuscript were complete. They left the protocol, endpoint, confirmatory results and primary inference unchanged and are labeled post hoc throughout.

Both controls reused the released prediction features, MPRA outcomes, variant metadata, region annotations, AVI scores and confirmatory utility results, inheriting eligibility, the emVar outcome, budgets of 384 and 1,000 variants per outcome cell and the SHA-256 tie-breaker. Absolute single-track quantiles were reconstructed row by row, with `magnitude_M` supplying the larger absolute quantile, `concordance_C` the smaller, and `magnitude_winning_modality` assigning them to DNase or H3K27ac. Across all 1,233,370 released rows, the reconstructed maximum, minimum and winning-modality identities had zero failures at an absolute tolerance of 1e-12. Combined-magnitude and AVI hit counts and precision exactly replayed all 16 corresponding confirmatory utility records before the new controls were evaluated.

The single-track control ranked each cell's candidates by DNase alone, H3K27ac alone, combined magnitude and AVI. The context control ranked the same candidate set five times using combined magnitude computed from A549, HepG2, K562, HCT116 or SK-N-SH predictions. The comparator was fixed as the unweighted mean precision of the four nonmatching rankings; no nonmatching context was selected on the basis of outcomes. Uncertainty used 2,000-replicate chromosome-stratified block bootstraps conditional on each fixed ranked selection at 1 Mb and 5 Mb, with seeds 21,261,008 and 25,261,008.

The released fields recover absolute rather than signed per-track quantiles, so these controls address prioritization precision and hit counts. Single-track directional accuracy lies outside their scope. The context control estimates performance relative to a fixed nonmatching mean across the five released contexts and does not causally decompose cell identity.

### Statistical analysis

We report cell-specific estimates and equal-weight confirmatory macro means. The frozen primary confidence interval used 2,000 chromosome-stratified 1 Mb block-bootstrap replicates (seed 20260913), with 5 Mb blocks (seed 20260918) for dependence sensitivity. Policies were reranked within each primary replicate; this remains the authoritative directional inference.

Fixed-budget intervals used 2,000 chromosome-stratified block replicates at 1 Mb and 5 Mb, conditional on the frozen selected lists. The same procedure estimated intervals for fixed top-10% rank-disagreement strata. These intervals characterize spatial sampling variability without policy reselection or optimization. Selective-risk curves at frozen coverage levels used analogous fixed-list intervals and remained descriptive because the primary reranking interval was already frozen.

Average precision and AUROC used each cell’s eligible fine-mapped universe; PIP metrics used observations with nonmissing PIP. Prespecified indicators were analyzed in separate cell-specific logistic regressions with 1 Mb cluster-robust standard errors. Benjamini–Hochberg FDR was controlled at 0.05 across 20 confirmatory-cell tests [6]. Interpretation rested on effects and confidence intervals rather than P values.

### Endogenous sensitivity evidence

We queried the versioned ADASTRA CL-TF BillCipher release distributed through IGVF [7]. Exact-cell allele-specific transcription-factor binding direction was `sign(effect_alt − effect_ref)`. Directional accuracy included variant–cell observations only when all significant transcription-factor observations agreed; ambiguous multi-factor observations were counted and excluded. This prespecified, mechanistically adjacent analysis could not redefine MPRA endpoints or policies.

### Reproducibility and integrity

Frozen SHA-256 hashes were `09bcfa0c09d6c39f6b98c521d73e772629395252f6faa4a35e4d3a8314d43f5f` for the protocol, `05df7123ab873d25b3d4b93d0e93dcfdeccad98b009c8fe55e73b39dff7f382c` for the pre-outcome evidence ledger and `f34b02fea83d59360ea4efadccfe9018c37a5910b8c033f1f69820bc2f03440e` for the confirmation-analysis script. The confirmation result inspection receipt records 25 of 25 checks passed. Manuscript analysis used Python 3.12.13 and NumPy 2.3.5 and passed 100 independent checks. The post hoc context-control verification receipt records 15 of 15 checks passed. All figures derive from source-data tables supplied with the repository.

## Data and code availability

The source MPRA and ADASTRA data are available from Siraj et al. [3] and IGVF [7], respectively. The public reproducibility repository provides the protocol, analysis code, derived data, figure source data and machine-readable Supplementary Tables S1–S9: https://github.com/LavindaLab/alphagenome-cell-matched-benchmark. The revised reproducibility package is versioned as `v1.1.0-biorxiv`; the earlier submitted version remains archived as release `v1.0.1-biorxiv`. AlphaGenome access is available from the provider [8].

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

### Figure 1 | Independent, outcome-blind evaluation of AlphaGenome Atlas

**A,** Frozen workflow from the Siraj MPRA through coordinate and allele quality control, Atlas retrieval, prediction freeze and held-out outcome analysis. **B,** Exact AlphaGenome ontology mapping and frozen inferential role for each MPRA cell line. Directional *n* is the number of noncoding, enhancer-active emVars available for directional analysis in confirmation blocks; A549 was prespecified as descriptive. **C,** Integrity checks. Atlas completeness refers to required DNase, H3K27ac and AVI scores. Prediction, result and manuscript-analysis inspection counts come from independent verification receipts. Gate outcomes were not used in confirmation, and no duplicate variant–cell rows were present.

### Figure 2 | Cell context improves AlphaGenome regulatory-variant prioritization

**A,** Equal-weight four-cell mean emVar precision for the fixed magnitude, concordance, AVI, PIP and deterministic pseudo-random policies at 96, 384, 1,000 and 5,000 assay slots per cell. Magnitude is the larger absolute exact-cell DNase/H3K27ac Atlas quantile score. **B,** Post hoc comparison of exact-cell macro precision with the unweighted mean of the four nonmatching prediction contexts at 384 and 1,000 slots. **C,** Post hoc exact-cell-minus-nonmatching-mean precision differences and 1 Mb block-bootstrap 95% confidence intervals by outcome cell. **D,** Additional emVars and additional 1 Mb blocks containing at least one emVar returned by exact-cell magnitude ranking relative to AVI across the four confirmatory cell-specific campaigns. **E,** Post hoc macro precision for AVI, DNase-only, H3K27ac-only and combined-magnitude ranking. Macro values are equal-weight means across HepG2, K562, HCT116 and SK-N-SH.

### Figure 3 | Disagreement is informative; concordance leaves ranking unimproved

**A,** EmVar prevalence among variants in the top 10% by cell-matched ranking but not AVI (“cell-matched only”) and by AVI but not cell-matched ranking (“AVI only”). Error bars are 1 Mb block-bootstrap 95% confidence intervals. **B,** Prevalence differences with 95% intervals; Macro is the equal-weight four-cell mean. **C,** Prespecified 10%-coverage sign-error contrast for concordance minus magnitude among enhancer-active emVars. Positive values favor magnitude. Intervals are from the frozen, chromosome-stratified 1 Mb bootstrap that reranked policies inside each replicate; this is the authoritative confirmatory analysis. **D,** Cluster-robust odds ratios for emVar status among DNase/H3K27ac-sign-agreeing versus sign-disagreeing variants. Agreement marks enrichment and leaves the primary directional ranking policy unimproved.

## References

1. Avsec Z, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature*. 2026;649:1206–1218. doi:10.1038/s41586-025-10014-0.
2. Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind technical report. 2026. Official landing page: https://deepmind.google/science/alphagenome/. Accessed 2026-09-14.
3. Siraj L, Castro RI, Dewey HB, et al. Functional dissection of complex trait variants at single-nucleotide resolution. *Nature*. 2026. doi:10.1038/s41586-026-10121-6.
4. Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. *Genome Research*. 2002;12:996–1006. doi:10.1101/gr.229102.
5. Frankish A, Diekhans M, Ferreira A-M, et al. GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Research*. 2019;47:D766–D773. doi:10.1093/nar/gky955.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289–300.
7. Abramov S, Boytsov A, Bykova D, et al. Landscape of allele-specific transcription factor binding in the human genome. *Nature Communications*. 2021;12:2751. doi:10.1038/s41467-021-23007-0.
8. Google DeepMind. AlphaGenome API README and terms summary. Accessed 2026-09-14. https://github.com/google-deepmind/alphagenome/blob/main/README.md.
