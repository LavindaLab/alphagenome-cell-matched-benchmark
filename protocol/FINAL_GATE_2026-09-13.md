# Final design gate: AlphaGenome Atlas external reliability study

**Date:** 2026-09-13
**Decision:** **MODIFY → conditional GO for a locked 10,000-SNV academic pilot**
**Full-corpus status:** not authorized until the pilot gate passes
**Implementation status:** no Atlas queries have been made

## 1. What the independent check learned

### Verified facts

1. **Generic abstention is not novel.** The original AlphaGenome paper selected score thresholds yielding 80% and 90% eQTL sign accuracy. At the 90% threshold it retained 41% of GTEx eQTLs, versus 19% for Borzoi, and it used the 80% threshold in GWAS credible sets.
2. **The Siraj MPRA corpus is usable at scale.** Official Nature Supplementary Table 4 contains 304,278 unique variants: exactly 281,327 canonical SNVs and 22,951 indels/noncanonical variants. It reports 92,560 active elements and 37,249 emVars. The article was published online on 2026-02-25.
3. **All five assay cell lines have exact Atlas ontology matches.** A549 (EFO:0001086), HepG2 (EFO:0001187), K562 (EFO:0002067), HCT116 (EFO:0002824) and SK-N-SH (EFO:0003072) occur in official AlphaGenome track metadata. DNase and H3K27ac tracks are present for every line. Other modalities differ by cell and cannot support a symmetric five-cell primary analysis.
4. **The primary data are rich enough without treating inactive constructs as equivalent to biological negatives.** Exact SNV counts are:

| Cell | Measured SNVs | Active SNVs | emVar SNVs | Active non-emVar SNVs |
|---|---:|---:|---:|---:|
| A549 | 114,451 | 8,719 | 2,574 | 6,145 |
| HCT116 | 173,166 | 28,495 | 6,782 | 21,713 |
| HepG2 | 280,409 | 55,329 | 15,184 | 40,145 |
| K562 | 280,097 | 55,782 | 18,200 | 37,582 |
| SK-N-SH | 279,999 | 53,319 | 16,114 | 37,205 |

5. **MPRA sign can be aligned to Atlas ALT-minus-REF sign.** The official MPRASuite documentation defines `Log2Skew = B_log2FC - A_log2FC` and labels A as reference and B as alternative. The Siraj variant identifier is `chromosome:position:allele1:allele2`, with allele 1 as reference in the inspected tables.
6. **Atlas bulk access has meaningful operational cost.** The official client issues one request per variant and defaults to ten workers. A full Siraj SNV run is therefore at most 281,327 requests, not five times that count if all ontology filters are included per request. No public fixed rate or latency SLA was found.
7. **The static AVI route is not a trivial download.** The current official terms page lists the permissive AVI SNV archive at 88.5 GB. It lists AVI feature importances at 283.9 GB and merged splicing scores at 20.6 GB in the noncommercial tier.
8. **Academic and business rights must remain separate.** Current terms allow academic/noncommercial API analysis and prohibit using outputs to train another ML model. The AVI SNV score archive is explicitly permissive for commercial and noncommercial use. Google Cloud documentation reviewed here covers the AlphaGenome base model, not Atlas; the Atlas launch describes commercial Atlas access as coming.

### Assumptions still requiring execution-time verification

- The supplied academic key has sufficient Atlas query capacity for the locked pilot.
- API scorer metadata returned at run time matches the official notebook snapshot.
- Siraj coordinates and alleles reconcile to GRCh38 at the required rate after reference validation.
- No newly released paper between this gate and posting occupies the narrow transportability claim.

## 2. What changed from the research brief

| Brief position | Corrected position |
|---|---|
| “When should Atlas abstain?” as the headline novelty | External **transportability** of score-based reliability plus incremental value of cross-modal agreement |
| Calibration language | Empirical, assay-conditional reliability; Atlas quantiles are not probabilities |
| MPRA as ground truth | MPRA as independent episomal evidence, with endogenous evidence needed for a central sensitivity result |
| Broad four-to-six-modality primary analysis | Symmetric five-cell primary analysis limited to exact-cell DNase and H3K27ac; other modalities secondary |
| 25,000-variant first pilot | 100 and 1,000 engineering checks, then a locked 10,000-SNV gate set |
| Possible raw-Atlas commercial routing through Cloud/customer environment | Not presently verified or allowed; Health AI is limited to permissive AVI, published aggregate findings and separately licensed inputs |
| Approximately 205,000 test SNVs | Exact all-corpus ceiling of 281,327 SNVs; the fine-mapped-only count is computed later from Table 1 rather than estimated |

## 3. Recommendation

Proceed with the frozen protocol in `PROTOCOL_v0.1.md`, using one paper with two linked aims:

1. **External reliability transport:** compare a magnitude-only Atlas policy with a fixed DNase–H3K27ac concordance policy at matched coverage across five exact cellular contexts.
2. **Experimental-budget utility:** replay fixed cell-specific selection policies at 96, 384, 1,000 and 5,000 assay slots, reporting hit yield and locus diversity.

Do not lead with AVI pathogenicity, generic benchmarking, generic abstention, “first independent validation,” or a learned uncertainty model.

## 4. Why this is the strongest near-term design

- It directly audits an already published scientific decision claim—high-confidence selection—under external assay and context shift.
- It exploits Atlas’s genome-wide precomputation while asking a question not answered by the AlphaGenome or Atlas papers.
- It has a negative-result path: stable non-transport or failure of cross-modal agreement to improve risk is scientifically informative.
- It uses all five exact cellular matches without pretending every modality is available everywhere.
- It preserves a clean academic–Health AI firewall and avoids an unsupported commercial Atlas assumption.
- It can stop after the pilot if the API, matching, or attainable precision makes a full 281,327-request run wasteful.

## 5. Next concrete action

Once `ALPHA_GENOME_API_KEY` is securely available, run the 100-variant metadata/schema check and persist the response manifest. If it passes, run 1,000 variants for rate/output sizing and then the locked 10,000-SNV gate set. Do not inspect the 90% confirmation set or download the 88.5 GB AVI archive before the gate decision.

The project is intentionally paused at this credential gate; this is an execution prerequisite, not a scientific access objection.
