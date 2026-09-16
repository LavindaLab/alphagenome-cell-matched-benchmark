# Held-out confirmation prequery freeze

**Frozen:** 2026-09-13, before confirmation MPRA outcome access  
**Query-set/code commit:** `5147add`  
**Pilot decision commit:** `a3d78f0`  
**Decision:** authorized to retrieve predictions; confirmation outcomes remain inaccessible

## Frozen scope

- 255,148 canonical SNVs occupy the prespecified 90% confirmation blocks.
- 246,674 remain after unique hg19→hg38 lift, reference validation, GENCODE 50 noncoding eligibility, primary-chromosome support and frozen prior-benchmark exclusions.
- Technical exclusions: 83 unmapped/nonunique lifts, 1,639 hg38 REF mismatches, 6,723 coding/CDS-or-essential-splice variants and 12 nonprimary-contig mappings.
- Seventeen exact Siraj candidates overlap an AlphaGenome CAGI5 interval or an explicit AlphaGenome/Atlas manuscript coordinate and are excluded before prediction retrieval.
- Required cell-matched scorer/track pairs are unchanged from the pilot: DNase track indices A549 107, HepG2 109, K562 118, HCT116 146, SK-N-SH 153; H3K27ac indices 328, 338, 392, 465 and 511 in the same order.
- One combined outcome-blind API request successfully returned `DNASE`, `CHIP_HISTONE` and `AVI_SCORE`; the confirmation run therefore uses one base request per variant rather than separate AVI requests.

## Query, time and storage projection

- Base requests: 246,674.
- Rate cap: 900 requests/minute with four workers and two bounded retries.
- Lower-bound retrieval time: 274.1 minutes (4.57 hours) before retries and service variance.
- Pilot-linear uncompressed estimates: 1.03 GB regulatory rows, 47 MB AVI rows and 224 MB derived features. A per-variant uncompressed cache would add about 1.51 GB; the confirmation implementation instead stores resumable deterministic gzip chunks, so expected committed summaries/features are much smaller while raw research-only chunks remain locally durable.
- No 88.5 GB static AVI archive is needed.

## Terms and version refresh

- Official SDK: `alphagenome` 0.9.0, Git commit `aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d`.
- Official Output Terms PDF inspected 2026-09-13: effective 2025-06-25, SHA-256 `a12935888c9c39e5568da86c5ec49fe651ad43ca3e37aa2184ad01b3e8320bc4`.
- The output terms permit free noncommercial use, prohibit commercial-organizational use except indirectly through scientific publication/open-source release, prohibit training a similar variant-effect model, require conspicuous downstream notice, and prohibit clinical use. This run is the academic/noncommercial lane. Raw outputs and derivatives are not Health AI commercial assets.
- Official README confirms Atlas predictions generally have higher query rates than on-demand model inference and describes the API as suitable below one million predictions; this run is 246,674 requests.

## Outcome firewall

The query manifest contains coordinates, alleles, block/split identity and measurement-availability metadata only. It contains no MPRA activity value, emVar call, skew, FDR, P value or label. The retrieval executable cannot read the Siraj workbooks. Raw cache records explicitly state that functional outcomes were not accessed.

Before confirmation outcomes are extracted, the final four-cell analysis executable, A549 downgrade, endogenous-evidence tier, scorer metadata receipt, retrieval inspection and all exact input hashes must be committed in a second freeze.
