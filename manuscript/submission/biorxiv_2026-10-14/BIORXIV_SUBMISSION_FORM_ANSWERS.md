# bioRxiv submission form answers

Prepared for the AlphaGenome manuscript on 14 October 2026.

## Core metadata

**Number of authors**

1

**Article category**

New Results

**Author type**

Select the portal's individual-author option.

**Subject area**

Genomics

**Title**

Cell-context matching improves AlphaGenome prioritization of experimentally active regulatory variants

Count: 10 words; 102 characters including spaces.

**Running title, if requested**

Cell-context AlphaGenome variant prioritization

Count: 4 words; 47 characters including spaces.

**Abstract - plain-text portal version**

AlphaGenome Atlas offers several predictions per variant, creating a choice for cell-specific experiments. We evaluated that choice against an independent massively parallel reporter assay (MPRA) covering 246,674 noncoding single-nucleotide variants in five human cell lines, with policies, endpoints and thresholds fixed before confirmatory outcomes were accessed. In a post hoc control that held the candidate set and scoring rule constant while changing only prediction context, exact-cell ranking reached 64.78% mean precision at 384 assay slots compared with 47.23% for the unweighted mean of four nonmatching contexts, a difference of 17.55 percentage points (1 Mb block-bootstrap 95% confidence interval, 15.35-19.67). At 1,000 slots the difference was 15.01 points (13.75-16.21). The benefit was positive in all four confirmatory cells, although it was small in SK-N-SH and some individual nonmatching contexts equaled or exceeded the exact context there. In the prespecified secondary comparison, exact-cell ranking reached 64.78% precision compared with 37.63% for the global AlphaGenome Variant Impact score (AVI), returning 417 additional expression-modulating variants across 367 additional megabase blocks at 384 slots and 973 across 674 blocks at 1,000. A post hoc single-track analysis found similar macro precision for H3K27ac alone and the combined DNase/H3K27ac score, while DNase alone was lower. The prespecified primary hypothesis failed: requiring DNase and H3K27ac effects to agree in sign made MPRA direction error 0.64 percentage points worse than magnitude ranking (-0.04 to 1.17). Across the tested cell lines, matching regulatory predictions to the assayed context improved experimental yield; cross-modal agreement annotated regulatory plausibility without improving directional ranking.

Count: 250 words. This plain-text copy replaces typographic range and minus symbols with ASCII characters for safe portal entry; its wording otherwise matches the manuscript abstract.

**Suggested keywords, if the portal requests them**

AlphaGenome; regulatory variant prioritization; MPRA; cell context; noncoding variants; variant effect prediction

## Author information

**Author and corresponding author**

Olga Lavinda

**Affiliation**

Department of Chemistry and Biochemistry, Stern College for Women, Yeshiva University, New York, NY 10016, USA

**Email**

olga.lavinda@yu.edu

**ORCID**

0000-0003-3577-5984

## Declarations and links

**Competing interests**

The author declares no competing interests.

**Funding**

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

If the funder field permits `No funder` or `Not applicable`, select it. Do not enter Yeshiva University as a funder unless institutional funds specifically supported this work.

**Ethics / human participants**

Not applicable. This computational study used publicly available variant-level assay data and did not analyze participant-level records or identifiable private information.

**Code and derived-data repository**

https://github.com/LavindaLab/alphagenome-cell-matched-benchmark

**Versioned release - verify this link after GitHub authentication and before final approval**

https://github.com/LavindaLab/alphagenome-cell-matched-benchmark/releases/tag/v1.1.0-biorxiv

**Source experimental data**

- Siraj MPRA article and associated supplementary data: https://doi.org/10.1038/s41586-026-10121-6
- IGVF ADASTRA reference file: https://data.igvf.org/reference-files/IGVFFI5943XCOS/
- AlphaGenome access and documentation: https://github.com/google-deepmind/alphagenome

**License**

CC BY 4.0

**Journal-submission status, if requested**

Answer from the live status on the day of bioRxiv submission. Concurrent journal submission is compatible with bioRxiv, but this worksheet does not infer whether a journal submission is currently active.

**Prior bioRxiv screening, only if directly requested**

The manuscript was previously submitted under BIORXIV/2026/751591 and BIORXIV/2026/752179 and was not posted. The present version adds a controlled exact-cell versus nonmatching-context analysis, DNase-only and H3K27ac-only controls, an updated Figure 2, and revised biological framing.

## Optional cover note

This manuscript reports a biological benchmarking study of cellular context in regulatory-variant prioritization. Using an independent MPRA across five human cell lines, it tests whether AlphaGenome regulatory predictions matched to the assayed cell improve recovery of experimentally active noncoding variants. A controlled analysis holds candidate sets and the ranking rule constant while changing only prediction context; exact-cell ranking improved precision relative to the fixed mean of four nonmatching contexts. The manuscript also reports the failed prespecified primary hypothesis and single-track controls. The public repository provides the protocol, analysis code, derived data, figure source data, machine-readable supplementary tables and verification receipts.

## Upload mapping

| Portal role | File | Suggested label |
| --- | --- | --- |
| Manuscript | `01_AlphaGenome_Cell_Context_bioRxiv_PREPRINT_2026-10-14.pdf` | Main manuscript |
| Supplemental file | `02_AlphaGenome_Supplementary_Methods_and_Results_2026-10-14.pdf` | Supplementary Methods and Results |
| Supplemental file | `03_AlphaGenome_bioRxiv_Supplementary_Data_2026-10-14.zip` | Supplementary Data 1: machine-readable Tables S1-S9 |

The manuscript PDF already contains all figures. Enter one manuscript file, zero separate image files and two supplemental files.

## Current portal facts and limits

Verified from bioRxiv's public submission help on 14 October 2026:

- Required stages include manuscript metadata, author information, file metadata, file upload, proof approval and confirmation.
- The online abstract must also appear in the manuscript file.
- One manuscript file is uploaded; a complete PDF with embedded figures requires no separate image files.
- Supplemental files are uploaded separately and should be labeled.
- Supplemental filenames should avoid symbols and punctuation other than hyphens, underscores and file extensions.
- The public help page does not state a numeric title, running-title or abstract limit.

The logged-in form may display counters or additional conditional questions that are not visible in public help. If it displays a limit below the counted values above, save the draft and stop before approval so the affected field can be shortened consistently in both the form and manuscript.

Official sources:

- Submission help: https://submit.biorxiv.org/help/submissionhelp.dtl
- Screening procedures: https://connect.biorxiv.org/news/2022/06/13/screening_procedures
- Funder fields: https://connect.biorxiv.org/news/2025/09/04/funder_information
- bioRxiv scope and licenses: https://www.biorxiv.org/about-biorxiv
