# Literature and claim lock

**Locked:** 2026-09-14 UTC, after the frozen results and before manuscript drafting
**Search scope:** exact-title and keyword searches through 2026-09-14 for AlphaGenome, AlphaGenome Atlas, AVI, MPRA, cell context, DNase, H3K27ac, functional-screen benchmarking and experimental prioritization
**Interpretation:** bounded novelty assessment, not proof that no unpublished or unindexed work exists

## Verified occupied ground

- AlphaGenome has already been evaluated on MPRA data. Avsec et al. used CAGI5 saturation MPRA, cell-matched DNase predictions, and learned cell-agnostic DNase and multimodal combinations. This manuscript therefore makes no “first AlphaGenome MPRA,” “first cell-matched,” “first multimodal,” or broad calibration claim.
- The official AlphaGenome Atlas preprint develops AVI from 18 global features, evaluates it on genome-wide and ten held-out saturation-genome-editing benchmarks, and demonstrates rare-disease and rare-variant prioritization. The official PDF inspected for this lock had SHA-256 `07011853872613cb7bd3ce8ca62cd2621ad50b0a16e6f5dec3e480b3d190af2c` (83 pages; manuscript creation date 2026-09-07). It contains no Siraj citation and no fixed-budget, exact-cell DNase/H3K27ac-versus-AVI MPRA analysis.
- Siraj et al. produced the independent five-cell MPRA evidence used here and analyzed activity, emVars, functional annotations and Enformer. Their paper contains no AlphaGenome analysis.
- Independent AlphaGenome work located by the lock includes personal-expression benchmarking, DNA-foundation-model benchmarks, A549/STARR-seq evaluation, distance-dependent calibration, MPRA fine-tuning, a one-locus RNU4-2/AVI comparison, a pan-cancer regulatory atlas, and disease-specific context filtering. These occupy broad external-validation, calibration and disease-prioritization territory, but not the present comparison.

## Narrow claim retained

The defensible contribution is an independent, outcome-blind comparison of two decision layers available in AlphaGenome Atlas: global AVI and exact-cell DNase/H3K27ac effects. In the Siraj five-cell MPRA corpus, the comparison asks which layer returns more experimentally measured regulatory variants under fixed assay budgets and whether disagreement marks failure of context-independent ranking.

This is narrower than model validation. The manuscript may claim:

1. a substantial, cross-cell experimental-yield advantage for the frozen exact-cell magnitude policy over AVI in this assay and candidate universe;
2. reproducibly lower MPRA emVar prevalence among AVI-only top-ranked variants than among matched-magnitude-only top-ranked variants at the prespecified 10% threshold; and
3. failure of the prespecified DNase/H3K27ac concordance policy to reduce MPRA direction error materially.

It may not claim universal superiority over AVI, causal or clinical validity, general probability calibration, in-vivo truth, or a “first validation.”

## Sources checked

- Avsec Z et al. *Nature* 649, 1206–1218 (2026). DOI: [10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0).
- Cheng J et al. *AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants*. Official launch preprint, released 2026-09-08: [PDF](https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf).
- Siraj L et al. *Nature* (2026). DOI: [10.1038/s41586-026-10121-6](https://doi.org/10.1038/s41586-026-10121-6).
- AlphaGenome official repository and terms summary, inspected 2026-09-14: [README](https://github.com/google-deepmind/alphagenome/blob/main/README.md).
- Nadal Vidal R. Preprints.org v1 (2026). DOI: [10.20944/preprints202607.1417.v1](https://doi.org/10.20944/preprints202607.1417.v1).
- Weykopf B et al. *PLOS Genetics* (2026). DOI: [10.1371/journal.pgen.1012222](https://doi.org/10.1371/journal.pgen.1012222).
- *AlphaGenome Enhances Personal Gene Expression Prediction but Retains Key Limitations*. bioRxiv v2 (2026): [10.1101/2025.08.05.668750](https://doi.org/10.1101/2025.08.05.668750).
- *Benchmarking DNA Foundation Models for Genomic and Genetic Tasks*. bioRxiv v2 (2025): [10.1101/2024.08.16.608288](https://doi.org/10.1101/2024.08.16.608288).
- *A pan-cancer regulatory atlas of 6,983 GWAS variants prioritizes recurrent regulatory annotations and candidate programs at cancer risk loci*. medRxiv v1 (2026): [10.64898/2026.05.16.26353369](https://doi.org/10.64898/2026.05.16.26353369).

## Posting-time action

Repeat the exact-title/keyword and forward-citation search on the day of public posting because Atlas was released only six days before this lock. Change positioning only if a new paper occupies the exact Siraj/exact-cell/AVI/fixed-budget comparison; do not alter the frozen analyses.
