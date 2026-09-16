# Substantive-value gate (pre-confirmation)

**Frozen:** 2026-09-13, after gate-pilot results and before confirmation outcomes
**Purpose:** decide whether results justify manuscript expansion; this is not a new inferential endpoint

Statistical significance alone cannot produce a GO. At least one candidate finding below must satisfy its quantitative screen and all five qualitative requirements. Otherwise the decision is **NO-GO**, the evaluation is preserved as a documented negative result, and work stops.

## Quantitative screens

### A. Experimental-prioritization utility

At both 384 and 1,000 slots per cell, matched DNase/H3K27ac magnitude must exceed AVI by at least 20% relative and 3 percentage points absolute in macro emVar precision across the four confirmatory cells, with the same direction in at least three of four cells. Report raw hit-count differences and locus diversity; do not substitute a P value for these thresholds.

### B. Concordance-policy effect

The four-cell primary sign-error difference must have absolute magnitude at least 2 percentage points, its 1 Mb block-bootstrap 95% interval must exclude zero, and the direction must agree in at least three of four cells. This may establish benefit or harm.

### C. Reliability-marker/ranking-policy distinction

DNase/H3K27ac sign agreement must associate with at least 1.5-fold emVar prevalence in at least three of four confirmatory cells, while the primary 1 Mb interval excludes a concordance benefit of 2 percentage points or more (lower bound greater than -0.02). This supports the narrower interpretation that agreement is useful for coarse reliability description but should not replace magnitude ranking.

## Five mandatory qualitative requirements

1. **Effect magnitude:** the effect crosses one of the practical screens above and is not driven by one cell or a tiny denominator.
2. **Novelty:** refreshed searches still distinguish the result from AlphaGenome CAGI5 MPRA benchmarking, Atlas AVI/SGE benchmarking, Weykopf A549 STARR-seq and Nadal Vidal cross-locus calibration.
3. **Experimental utility:** a scientist can make a concrete assay-allocation or abstention/ranking choice differently because of the result.
4. **Interpretability:** the direction is mechanistically legible, survives 5 Mb dependence sensitivity where applicable, and is not contradicted without explanation by exact-cell ADASTRA allele-specific binding.
5. **Change to use:** the conclusion specifies when to prefer matched molecular scores, AVI, concordance or no automated prioritization; it cannot be merely “AlphaGenome is statistically associated with MPRA.”

## Decision rule

- **GO:** at least one quantitative screen plus all five qualitative requirements. Manuscript expansion may then begin, with the gate result and limitations stated explicitly.
- **NO-GO:** no candidate passes. Preserve hashes, code, tables, intervals, novelty record and a concise negative-result report; perform no manuscript expansion or optional analyses.

The gate is deliberately stricter than statistical significance and was written without confirmation outcome access. It does not alter the frozen primary population, policies, coverage, outcome or estimand.
