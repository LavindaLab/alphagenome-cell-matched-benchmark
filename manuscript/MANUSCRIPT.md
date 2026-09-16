# Cell-matched AlphaGenome predictions improve experimental variant prioritization over a global impact score

**Olga Lavinda**

Department of Chemistry and Biochemistry, Stern College for Women, Yeshiva University, New York, NY 10016, USA

Correspondence: olga.lavinda@yu.edu · ORCID: 0000-0003-3577-5984

## Abstract

AlphaGenome Atlas provides several predictions for each variant, leaving scientists to choose between a global impact score and molecular outputs matched to the planned experiment. In a prespecified, outcome-blind analysis of an independent massively parallel reporter assay (MPRA), cell-matched DNase/H3K27ac ranking achieved 64.78% mean precision at 384 assay slots per confirmatory cell versus 37.63% for the global AlphaGenome Variant Impact score (AVI), a difference of 27.15 percentage points (1 Mb block-bootstrap 95% confidence interval, 22.93–31.22). At 1,000 slots, the difference was 24.32 points (21.87–26.82). The analysis included 246,674 noncoding single-nucleotide variants across five human cell lines; four were confirmatory and A549 was prespecified as descriptive. Cell-matched ranking returned 417 and 973 additional expression-modulating variants at the two budgets and performed better in every confirmatory context. Among variants ranked in the top 10% by only one score, equal-weight mean hit prevalence was 14.71% for cell-matched-only variants and 6.33% for AVI-only variants (difference 8.38 points, 7.84–8.91); the equal-weight mean of the four cell-specific prevalence ratios was 2.41 (2.27–2.57). The study’s prespecified primary hypothesis was different and failed: requiring DNase and H3K27ac effects to agree in sign made MPRA direction error 0.64 points worse than magnitude ranking (−0.04 to 1.17). In this assay and candidate universe, matching the prediction layer to cellular context materially improved experimental return; cross-modal agreement annotated regulatory plausibility but did not improve directional ranking.

## Introduction

AlphaGenome Atlas creates a practical choice for experimental design: when several predictions are available for the same variant, which should determine what is tested? Atlas provides molecular effects for specific cell types and assays, alongside AVI, a global score that combines AlphaGenome outputs with conservation and protein-related features [1,2]. The Atlas interface and launch materials foreground AVI as the single-number route for rapid variant ranking [2], making it the lowest-friction starting point for candidate selection. Genome-wide precomputation removes much of the inference burden from variant screening, but it does not resolve whether a broad score or a biologically matched output is the better guide to a context-specific experiment.

AVI supports broad prioritization across coding and noncoding mechanisms. Exact-cell DNase and H3K27ac predictions instead represent regulatory processes closer to the evidence produced by an MPRA. AlphaGenome’s original study benchmarked MPRA effects, including matched DNase predictions and learned combinations [1]; the Atlas study benchmarked AVI broadly and against saturation genome-editing screens [2]. Neither tested the zero-shot choice between AVI and exact-cell Atlas outputs on the large five-cell MPRA reported by Siraj and colleagues [3].

We used that independent MPRA for a frozen, outcome-blind evaluation. The primary hypothesis asked whether requiring DNase and H3K27ac effects to agree in sign would reduce directional error at fixed coverage. A prespecified secondary aim compared experimental yield under frozen ranking policies and assay budgets. The primary hypothesis failed. Cell-matched ranking nevertheless produced substantially more MPRA hits than AVI, and rank disagreement identified variants for which AVI was a poor guide to this context-specific assay. We center the paper on this actionable prespecified secondary comparison and retain the failed primary as a central result. Together, the results separate two decisions: match the prediction layer to the experiment, but do not treat cross-modal agreement as evidence of more reliable reporter-assay direction.

## Results

### An outcome-blind, independent evaluation separated engineering from confirmation

We used the 304,278-variant MPRA of Siraj et al., which tested fine-mapped trait and expression quantitative trait locus variants and matched controls in K562, HepG2 and SK-N-SH, with subsets in A549 or HCT116 [3]. We retained the source authors’ cell-specific definition of an expression-modulating variant (emVar): an active element with significant allelic skew at false-discovery rate below 0.1.

Policy definitions, contexts, scores, coverage thresholds, assay budgets and outcome definitions were frozen before confirmation-outcome access (Fig. 1). Genomic blocks were assigned deterministically to a 10% gate set or 90% confirmation set. Gate outcomes informed only the decision to run the held-out analysis and were permanently excluded from confirmation. A prespecified event-count rule retained HepG2, K562, HCT116 and SK-N-SH for equal-weight confirmatory summaries; A549 remained descriptive.

Source hg19 variants were uniquely mapped to hg38 and checked against the reference allele. We excluded invalid or nonprimary mappings, reference mismatches, protein-coding CDS or essential-splice variants under GENCODE 50, and 17 exact overlaps with identified AlphaGenome or Atlas benchmark/example coordinates. The resulting 246,674 unique noncoding SNVs had complete predictions across all five exact-cell ontologies (1,233,370 variant–cell rows). Independent inspections passed every retrieval, frozen-result and manuscript-analysis check; no duplicate row or gate/confirmation crossover was detected.

### Cell-matched DNase/H3K27ac magnitude substantially outperformed AVI at fixed budgets

The cell-matched magnitude policy ranked variants by the larger absolute quantile score from the exact-cell DNase and H3K27ac predictions. AVI ranked the same eligible fine-mapped, noncoding candidates by the Atlas global-score quantile. Neither policy used MPRA outcomes for construction or threshold selection.

At 384 selected variants per cell, cell-matched magnitude achieved a four-cell macro precision of 0.6478 versus 0.3763 for AVI, an absolute difference of 0.2715 (1 Mb block-bootstrap 95% CI, 0.2293–0.3122; Fig. 2A,B). At 1,000 variants, precision was 0.5625 versus 0.3193, a difference of 0.2432 (0.2187–0.2682; Fig. 2A,C). The advantage was positive in every confirmatory context (Supplementary Table S1).

Cell-matched ranking remained better than AVI, the frozen PIP baseline and random selection at all four budgets (96, 384, 1,000 and 5,000; Supplementary Fig. S1 and Supplementary Table S1). The macro precision difference over AVI declined from 0.4401 at 96 slots (0.3550–0.5202) to 0.1253 at 5,000 (0.1170–0.1336). AVI precision rose between 96 and 384 slots because ranks 97–384 had higher aggregate emVar prevalence than its top 96; precision at a larger cutoff is not constrained to decline. Five-megabase dependence blocks gave nearly identical intervals (Supplementary Fig. S2). Average precision and AUROC also favored cell-matched ranking in each cell (Supplementary Table S2).

### Rank disagreement identifies a reproducible loss of experimental yield

At the prespecified 10% threshold, we classified variants as highly ranked by both policies, cell-matched only, AVI only or neither. This descriptive comparison of frozen scores was not an optimized abstention rule. Across confirmatory contexts, equal-weight mean emVar prevalence was 0.1471 among cell-matched-only variants and 0.0633 among AVI-only variants, a difference of 0.0838 (1 Mb block-bootstrap 95% CI, 0.0784–0.0891). The equal-weight mean of the four cell-specific prevalence ratios was 2.41 (2.27–2.57; Fig. 3A,B); it is distinct from the ratio of the two macro prevalences.

The separation was present in every confirmatory context; cluster-robust odds ratios for cell-matched-only versus AVI-only membership ranged from 1.61 to 3.93 after Benjamini–Hochberg control across the 20 prespecified indicator tests (Supplementary Fig. S3 and Supplementary Tables S3, S5).

DNase–H3K27ac sign disagreement, a larger cross-modal quantile gap and exact-cell versus other-cell direction disagreement were also associated with lower emVar odds in all four contexts. In this dataset, disagreement therefore functions as a warning annotation, not a calibrated uncertainty probability or portable decision threshold.

### The prespecified concordance hypothesis failed

The primary endpoint was MPRA allelic-direction error at 10% coverage among enhancer-active emVars. Magnitude predicted direction from the modality with the larger effect. Concordance required nonzero DNase and H3K27ac effects with the same sign and ranked eligible observations by the smaller absolute quantile score.

The hypothesis failed. Four-cell macro sign error was 0.64 percentage points worse under concordance than magnitude ranking (prespecified contrast, +0.00637; reranking 1 Mb block-bootstrap 95% CI, −0.00040 to 0.01168; Fig. 3C). The 5 Mb interval was −0.00004 to 0.01115. Both intervals included zero, and the point estimate was well below the prespecified 2-point threshold for practical importance. Cell-specific estimates had the same direction, while fixed-list risk–coverage differences were small, coverage-dependent and changed sign (Supplementary Table S4). No MPRA analysis supported a general concordance advantage.

Agreement marked a different property: emVar prevalence was 1.39–1.76 times higher among sign-agreeing than sign-disagreeing variants, with corresponding logistic-regression odds ratios of 1.41–1.82 (Fig. 3D). In exact-cell ADASTRA allele-specific transcription-factor binding, concordance modestly improved direction accuracy in four lines and was nearly neutral in SK-N-SH. This mechanistically adjacent evidence concerns endogenous binding, not reporter transcription, and does not alter the negative MPRA endpoint. Agreement can annotate regulatory plausibility without improving directional ranking.

### Matched ranking increased experimentally confirmed discoveries and locus diversity

At 384 slots per cell, cell-matched ranking returned 995 emVars versus 578 for AVI: 417 additional hits across 367 additional hit-containing 1 Mb blocks (Fig. 2D). At 1,000 slots, it returned 2,250 versus 1,277 hits, a gain of 973 across 674 additional hit-containing blocks. The improvement was therefore not confined to repeatedly selecting nearby variants at a few loci.

For a noncoding reporter experiment in these contexts, exact-cell accessibility and active-histone effects produced more confirmed regulatory variants per synthesis slot than AVI. This is not a universal comparison of score quality: AVI addresses broader coding and noncoding impact and remains appropriate when context is unknown, a relevant track is unavailable or the intended evidence is not regulatory activity.

## Discussion

AlphaGenome Atlas makes generic and context-specific predictions equally accessible, but they are not interchangeable experimental guides. In this independent MPRA, ranking noncoding variants by exact-cell DNase/H3K27ac effect improved hit precision over AVI by 24–27 percentage points at two realistic assay budgets. It returned hundreds of additional hits across hundreds of additional loci, with consistent gains across liver-, blood-, colon- and neural-derived cell lines and two genomic block sizes. For a matched regulatory assay, scientists should rank the relevant molecular outputs rather than defaulting to AVI.

The prespecified primary hypothesis failed. DNase and H3K27ac agreement enriched for emVars but did not reduce reporter allelic-direction error relative to magnitude ranking. Agreement records whether two regulatory layers point in the same direction; it does not establish that their shared sign will carry through a 200-bp episomal reporter. Requiring agreement can also exclude strong single-modality effects. The small ADASTRA sensitivity analysis was compatible with modest directional information for endogenous binding, but its different endpoint cannot resolve the negative MPRA result. In this study, agreement is an annotation, not a better ranking policy.

Rank disagreement makes the distinction operational. Variants selected only by matched regulatory magnitude had more than twice the emVar prevalence of variants selected only by AVI at the same 10% threshold. An AVI-high, matched-low variant should not therefore be assumed to be a high-yield reporter candidate in these settings. AVI may reflect conservation, coding, splicing or cross-context effects that the MPRA does not measure; disagreement indicates misalignment with this assay, not that AVI is biologically wrong.

The inference is bounded by the assay and candidate set. MPRA measurements are episomal and assay-conditional: they provide independent experimental evidence, not universal in-vivo ground truth. Candidates were enriched for fine-mapped human variants and controls rather than sampled uniformly from the genome. Exact matching required cell lines with both DNase and H3K27ac tracks, so transfer to primary tissues, perturbations or uncovered contexts is unknown. Because AVI and the molecular scores come from the same model resource, this study compares decision layers rather than independent models. Fixed-budget intervals condition on frozen selected lists; the primary directional interval reranked within bootstrap replicates and remains authoritative. Novelty is bounded by the literature search through 2026-09-14.

The reported effect belongs to the prespecified maximum-magnitude combination of DNase and H3K27ac. Attribution of the gain to either component track requires a separately declared analysis.

Prediction layers should be selected for the evidence an experiment will produce. For cell-specific noncoding reporter assays, context matching materially changed experimental yield. A plausible cross-modal concordance requirement did not improve directional reliability.

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

### Statistical analysis

We report cell-specific estimates and equal-weight confirmatory macro means. The frozen primary confidence interval used 2,000 chromosome-stratified 1 Mb block-bootstrap replicates (seed 20260913), with 5 Mb blocks (seed 20260918) for dependence sensitivity. Policies were reranked within each primary replicate; this remains the authoritative directional inference.

Fixed-budget intervals used 2,000 chromosome-stratified block replicates at 1 Mb and 5 Mb, conditional on the frozen selected lists. The same procedure estimated intervals for fixed top-10% rank-disagreement strata. These intervals characterize spatial sampling variability without policy reselection or optimization. Selective-risk curves at frozen coverage levels used analogous fixed-list intervals and remained descriptive because the primary reranking interval was already frozen.

Average precision and AUROC used each cell’s eligible fine-mapped universe; PIP metrics used observations with nonmissing PIP. Prespecified indicators were analyzed in separate cell-specific logistic regressions with 1 Mb cluster-robust standard errors. Benjamini–Hochberg FDR was controlled at 0.05 across 20 confirmatory-cell tests [6]. Interpretation rested on effects and confidence intervals rather than P values.

### Endogenous sensitivity evidence

We queried the versioned ADASTRA CL-TF BillCipher release distributed through IGVF [7]. Exact-cell allele-specific transcription-factor binding direction was `sign(effect_alt − effect_ref)`. Directional accuracy included variant–cell observations only when all significant transcription-factor observations agreed; ambiguous multi-factor observations were counted and excluded. This prespecified, mechanistically adjacent analysis could not redefine MPRA endpoints or policies.

### Reproducibility and integrity

Frozen SHA-256 hashes were `09bcfa0c09d6c39f6b98c521d73e772629395252f6faa4a35e4d3a8314d43f5f` for the protocol, `05df7123ab873d25b3d4b93d0e93dcfdeccad98b009c8fe55e73b39dff7f382c` for the pre-outcome evidence ledger and `f34b02fea83d59360ea4efadccfe9018c37a5910b8c033f1f69820bc2f03440e` for the confirmation-analysis script. Confirmation result commit `0422dbf7678f0c84fda9d6a9336fc7b1b10ad8e2` passed 25 independent checks. Manuscript analysis used Python 3.12.13 and NumPy 2.3.5 and passed 100 independent checks. All figures derive from committed tables supplied as source data.

## Data and code availability

The source MPRA and ADASTRA data are available from Siraj et al. [3] and IGVF [7], respectively. The public reproducibility repository provides the protocol, analysis code, derived data, figure source data and machine-readable Supplementary Tables S1–S7: https://github.com/healthai-hq/alphagenome-cell-matched-benchmark. The exact preprint release is archived at https://github.com/healthai-hq/alphagenome-cell-matched-benchmark/releases/tag/v1.0.0-biorxiv. AlphaGenome access is available from the provider [8].

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

### Figure 2 | Cell-matched regulatory scores improve experimental yield

**A,** Equal-weight four-cell mean emVar precision for the frozen magnitude, concordance, AVI, PIP and deterministic pseudo-random policies at 96, 384, 1,000 and 5,000 assay slots per cell. Magnitude is the larger absolute exact-cell DNase/H3K27ac Atlas quantile score. Concordance is a distinct policy from magnitude; their prespecified directional comparison is shown in Figure 3C. **B,C,** Magnitude-minus-AVI precision differences and 95% confidence intervals at 384 and 1,000 slots, estimated by 2,000 chromosome-stratified 1 Mb block-bootstrap replicates conditional on each frozen selected list. Macro is the equal-weight mean across HepG2, K562, HCT116 and SK-N-SH. Dashed lines mark no difference. **D,** Additional emVars and additional 1 Mb blocks containing at least one emVar returned by magnitude ranking across the four confirmatory cell-specific campaigns.

### Figure 3 | Disagreement is informative, but concordance is not a better ranking rule

**A,** EmVar prevalence among variants in the top 10% by cell-matched ranking but not AVI (“cell-matched only”) and by AVI but not cell-matched ranking (“AVI only”). Error bars are 1 Mb block-bootstrap 95% confidence intervals. **B,** Prevalence differences with 95% intervals; Macro is the equal-weight four-cell mean. **C,** Prespecified 10%-coverage sign-error contrast for concordance minus magnitude among enhancer-active emVars. Positive values favor magnitude. Intervals are from the frozen, chromosome-stratified 1 Mb bootstrap that reranked policies inside each replicate; this is the authoritative confirmatory analysis. **D,** Cluster-robust odds ratios for emVar status among DNase/H3K27ac-sign-agreeing versus sign-disagreeing variants. Agreement marks enrichment but did not improve the primary directional ranking policy.

## References

1. Avsec Z, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature*. 2026;649:1206–1218. doi:10.1038/s41586-025-10014-0.
2. Cheng J, Taylor KR, Nicolaisen L, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. Google DeepMind technical report. 2026. Official landing page: https://deepmind.google/science/alphagenome/. Accessed 2026-09-14.
3. Siraj L, Castro RI, Dewey HB, et al. Functional dissection of complex trait variants at single-nucleotide resolution. *Nature*. 2026. doi:10.1038/s41586-026-10121-6.
4. Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. *Genome Research*. 2002;12:996–1006. doi:10.1101/gr.229102.
5. Frankish A, Diekhans M, Ferreira A-M, et al. GENCODE reference annotation for the human and mouse genomes. *Nucleic Acids Research*. 2019;47:D766–D773. doi:10.1093/nar/gky955.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289–300.
7. Abramov S, Boytsov A, Bykova D, et al. Landscape of allele-specific transcription factor binding in the human genome. *Nature Communications*. 2021;12:2751. doi:10.1038/s41467-021-23007-0.
8. Google DeepMind. AlphaGenome API README and terms summary. Accessed 2026-09-14. https://github.com/google-deepmind/alphagenome/blob/main/README.md.
