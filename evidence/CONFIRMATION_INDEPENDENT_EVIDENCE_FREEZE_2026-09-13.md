# Confirmation independent-evidence freeze

**Frozen:** 2026-09-13, before held-out MPRA outcome access
**Purpose:** provide independent biological direction evidence without using Siraj MPRA labels

## Tier 1 attempt: IGVF Variant-EFFECTS

The official IGVF Catalog production API returned 737 records on 2026-09-13. All were prime-editing records and none were SNVs, so this source cannot independently validate the frozen SNV study. The complete response and receipt are preserved rather than silently substituting it after outcomes are known.

- API response: `igvf_variant_effects_2026-09-13.jsonl.gz`
- Receipt: `igvf_variant_effects_2026-09-13_receipt.json`
- Response SHA-256: `22680b11500e769603257463dca514a720290fcc7d3927c59f8474ce4fd1bbba`

## Prespecified Tier 2 fallback: exact-cell allele-specific TF binding

Official IGVF reference file `IGVFFI5943XCOS` distributes the ADASTRA CL-TF BillCipher July 2022 release. The hg38 archive was intersected by exact allele and exact cell name with the frozen confirmation manifest.

- 378 exact-cell TF files inspected.
- 31,266 archive rows scanned across A549, HepG2, K562, HCT-116 and SK-N-SH.
- 3,340 overlapping rows, all meeting the archived 5% FDR indicator under the frozen extraction rule.
- 2,329 unique variant-cell pairs.
- Direction is `sign(es_mean_alt - es_mean_ref)`; agreement across significant TFs is required for a variant-cell consensus.
- Output SHA-256: `e199d3eddda863f3805fc1412f3fb7640a5ef0a797b5d747c6cb0ea63a8ad856`.
- Nested ADASTRA archive SHA-256: `d8825eaf985baa4d8cea7b5215a4cfea988065abcc12022b88f911a7d763f648`.

The sensitivity is secondary and mechanistically limited: TF-binding direction is not equivalent to enhancer activity or MPRA allelic skew. It is used only to assess whether an apparent ranking/direction result has independent exact-cell support or unexplained contradiction. It cannot tune any rule.

## Firewall

Neither extraction opened Siraj Supplementary Table 4. Both receipts record that confirmation MPRA outcomes were not accessed and that scientific rules were not mutated. The original `EVIDENCE_LEDGER.md` remains unmodified as the preserved pre-outcome ledger; this file is a dated confirmation supplement.
