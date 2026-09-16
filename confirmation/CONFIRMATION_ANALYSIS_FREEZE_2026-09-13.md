# Held-out confirmation analysis freeze

**Frozen:** 2026-09-13, before confirmation MPRA outcome access
**Outcome state:** inaccessible; neither Table 4 functional labels nor confirmation outcome rows have been opened
**Scientific rules mutated:** no

## Confirmatory analysis

- The primary population, outcome, policies, directions, coverage and estimand remain those in `PROTOCOL_v0.1.md`.
- The prespecified pilot count rule is active: HepG2, K562, HCT116 and SK-N-SH form the equal-weight confirmatory macro-average; A549 is descriptive only.
- The primary contrast is `sign-error(concordant) - sign-error(magnitude)` at 10% coverage among measured, Atlas-complete, noncoding, enhancer-active emVars. Positive values favor magnitude ranking.
- Confidence intervals use 2,000 chromosome-stratified 1 Mb block-bootstrap replicates with seed `20260913`. The 5 Mb dependence sensitivity uses seed `20260918`.
- Risk-coverage estimates remain fixed at 1%, 5%, 10%, 25% and 50%.

## Experimental-prioritization analysis

- The population is the measured, Atlas-complete, noncoding, fine-mapped candidate universe in each cell.
- Policies are deterministic chromosome/PIP-availability-stratified random, PIP, AVI, matched DNase/H3K27ac magnitude, matched DNase/H3K27ac concordance, and concordance with a maximum of two selections per 1 Mb block.
- Budgets are 96, 384, 1,000 and 5,000 variants per cell.
- Outcomes are raw emVar hits, precision, cumulative absolute MPRA skew, independent 1 Mb blocks with a hit, missed emVars, and the promoter-proximal composition of selected variants and hits.
- Promoter-proximal is frozen as within 2,000 bp of a GENCODE 50 protein-coding transcript TSS. It is descriptive and cannot change rankings or primary endpoints.

## Discordance and independent-evidence analyses

- DNase/H3K27ac sign agreement is compared with disagreement using emVar prevalence and magnitude-policy directional error among enhancer-active emVars.
- The independent sensitivity is exact-cell allele-specific TF binding from the ADASTRA CL-TF BillCipher release distributed by IGVF. Direction is `sign(es_mean_alt - es_mean_ref)`. Variant-cell evidence is retained only when all significant TF observations agree in direction; ambiguous multi-TF observations are counted and excluded from directional accuracy.
- This evidence is mechanistically adjacent rather than an MPRA substitute. It may support or challenge interpretation but cannot redefine thresholds, policies, populations or endpoints.

## Substantive-value decision

`SUBSTANTIVE_VALUE_GATE_v1.md` is frozen before outcome access. A quantitative screen is necessary but insufficient: a GO also requires meaningful magnitude, novelty after a post-result literature refresh, concrete experimental utility, mechanistic interpretability including the 5 Mb and ADASTRA checks, and a clear change to how a scientist should use AlphaGenome. Failure is a documented NO-GO and stops the project before manuscript expansion.

## Frozen outcome-blind inputs

| Artifact | Rows/role | SHA-256 |
|---|---:|---|
| `variants_hg38.tsv.gz` | 246,674 confirmation SNVs | `10761f6ef1f4c46f31175fbdcdaf3f119943a7c6fb9a951864d6c4766ff986a1` |
| `variant_metadata.tsv.gz` | PIP/category baseline metadata | `f3f6aa76e2a0468c9c455b3d930b16db457f42f19b4d927a94b16a7c9aae8681` |
| `tss_distance.tsv.gz` | GENCODE 50 TSS annotation | `62f1d96b04d30a58578120b46ac17cf4ef844c13c2937d101a1591620b334125` |
| `adastra_exact_cell_confirmation.tsv.gz` | 3,340 TF-by-variant-cell overlap rows; 2,329 unique variant-cell pairs | `e199d3eddda863f3805fc1412f3fb7640a5ef0a797b5d747c6cb0ea63a8ad856` |
| `igvf_variant_effects_2026-09-13.jsonl.gz` | 737 Tier-1 records; no SNVs, therefore unsuitable | `22680b11500e769603257463dca514a720290fcc7d3927c59f8474ce4fd1bbba` |
| `analyze_confirmation.py` | held-out analysis executable | `f34b02fea83d59360ea4efadccfe9018c37a5910b8c033f1f69820bc2f03440e` |

Prediction-feature, AVI, cache-tree and inspection hashes will be appended through the separate outcome-access authorization only after retrieval is complete and independently inspected. That later receipt cannot alter this analysis or the substantive-value gate.
