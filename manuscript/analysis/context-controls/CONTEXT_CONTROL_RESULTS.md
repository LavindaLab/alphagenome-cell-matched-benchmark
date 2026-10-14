# Post hoc AlphaGenome context controls

## Scope

These analyses are post hoc diagnostics requested on 8 October 2026. They leave the protocol, confirmatory endpoint, original result record, and manuscript unchanged. The confirmatory outcome cells are HepG2, K562, HCT116, and SK-N-SH; A549 is used only as a nonmatching prediction context.

## Methods

The analyses used the released prediction features, MPRA outcomes, variant metadata, region annotations, AVI scores, and frozen utility results. Eligibility, the emVar outcome, budgets of 384 and 1,000 variants per outcome cell, and the SHA-256 tie breaker were inherited from the original analysis. Eligible populations were 176,396 variants for HepG2, 176,198 for K562, 119,946 for HCT116, and 176,100 for SK-N-SH.

Absolute single-track quantiles were reconstructed row by row. `magnitude_M` supplied the larger absolute quantile, `concordance_C` supplied the smaller, and `magnitude_winning_modality` assigned those values to DNase or H3K27ac. Across all 1,233,370 released prediction-feature rows, the reconstructed maximum, minimum, and winning-modality identities had zero failures at an absolute tolerance of 1e-12. A row-level example is retained in `context_control_summary.json`.

For each confirmatory outcome cell, DNase-only, H3K27ac-only, combined-magnitude, and AVI rankings used the same eligible candidates. The exact-versus-wrong-cell control ranked that same candidate set five times using combined magnitude from A549, HepG2, K562, HCT116, or SK-N-SH. The fixed wrong-cell comparator was the unweighted mean precision of the four nonmatching rankings. No wrong context was selected based on outcomes.

Uncertainty used the repository's established 2,000-replicate chromosome-stratified block bootstrap, conditional on each fixed ranked selection, at 1 Mb and 5 Mb. Seeds were 21,261,008 and 25,261,008, respectively. Combined-magnitude and AVI hit counts and precision exactly replayed all 16 corresponding frozen utility records before the new controls were evaluated.

## Results

### Single-track control

Macro precision is the unweighted mean over the four confirmatory outcome cells. Hit counts are totals over those four equally sized selections.

| Budget | DNase-only | H3K27ac-only | Combined magnitude | AVI |
|---:|---:|---:|---:|---:|
| 384 | 0.5716 (878 hits) | 0.6673 (1,025 hits) | 0.6478 (995 hits) | 0.3763 (578 hits) |
| 1,000 | 0.5155 (2,062 hits) | 0.5590 (2,236 hits) | 0.5625 (2,250 hits) | 0.3193 (1,277 hits) |

H3K27ac-only captured the combined score's macro performance. At budget 384 it was 0.0195 higher than combined magnitude (1 Mb 95% interval, -0.0041 to 0.0442; 5 Mb, -0.0034 to 0.0413) and produced 30 more hits. At budget 1,000 it was 0.0035 lower (1 Mb, -0.0157 to 0.0082; 5 Mb, -0.0146 to 0.0087) and produced 14 fewer hits. Thus H3K27ac-only exceeded the combined-minus-AVI point gain at budget 384 and retained 98.6% of that gain at budget 1,000. The combined maximum showed no clear macro improvement over H3K27ac-only in either block analysis.

DNase-only retained useful prioritization signal but was lower than combined magnitude by 0.0762 at budget 384 (1 Mb, -0.0963 to -0.0567; 5 Mb, -0.0949 to -0.0583) and by 0.0470 at budget 1,000 (1 Mb, -0.0586 to -0.0358; 5 Mb, -0.0587 to -0.0362). The cell-level pattern was heterogeneous: H3K27ac-only was particularly favorable in SK-N-SH, whereas the combined score was higher in K562 and HCT116 at budget 1,000.

### Exact-cell versus all-wrong-cell control

| Budget | Exact-cell macro precision | All-wrong-cell mean precision | Exact minus all-wrong mean | 1 Mb 95% interval | 5 Mb 95% interval |
|---:|---:|---:|---:|---:|---:|
| 384 | 0.6478 | 0.4723 | 0.1755 | 0.1535 to 0.1967 | 0.1536 to 0.1979 |
| 1,000 | 0.5625 | 0.4124 | 0.1501 | 0.1375 to 0.1621 | 0.1373 to 0.1632 |

The fixed exact-minus-all-wrong contrast was positive in every outcome cell and at both budgets. At budget 384, the 1 Mb cell-level differences were 0.0892 for HepG2, 0.3737 for K562, 0.1895 for HCT116, and 0.0495 for SK-N-SH; their respective 95% intervals were 0.0355 to 0.1409, 0.3293 to 0.4190, 0.1525 to 0.2276, and 0.0067 to 0.0908. At budget 1,000 the differences were 0.0988, 0.3215, 0.1518, and 0.0283, with intervals of 0.0681 to 0.1303, 0.2949 to 0.3490, 0.1289 to 0.1748, and 0.0003 to 0.0521. The 5 Mb sensitivity analysis preserved these conclusions.

Exact-cell ranking was the highest of the five prediction contexts for HepG2, K562, and HCT116. For SK-N-SH, some individual wrong contexts matched or exceeded the exact context even though the prespecified fixed all-wrong mean was lower. The control therefore supports an average benefit of exact cell matching across the four outcome cells, with a smaller and less individually dominant benefit in SK-N-SH.

## Limitations

These diagnostics were specified after the confirmatory results and should remain labeled post hoc. The released fields recover absolute single-track quantiles but not the signed per-track quantiles, so the control addresses prioritization precision and hit counts rather than single-track directional accuracy. The block bootstrap conditions on each fixed ranked selection. The wrong-cell comparison covers the five released cellular contexts and estimates performance relative to their fixed nonmatching mean; it is not a causal decomposition of cell identity or a claim that the exact context beats every individual mismatch.

## Conditional manuscript-framing guide

- Exact-cell context: the fixed control supports stating that exact-cell combined-magnitude ranking improved average emVar precision relative to the unweighted mean of all four nonmatching contexts at both budgets. Keep the scope to these four outcome cells and the five tested prediction contexts, and avoid universal claims about every individual mismatch.
- Two-track maximum: H3K27ac-only accounts for essentially all of the combined score's macro prioritization gain over AVI at both budgets. The post hoc results do not support attributing the gain specifically to taking the maximum across DNase and H3K27ac. DNase-only contributes signal but does not reproduce the combined macro performance.
- Decision rule: the observed exact-cell benefit is consistently positive under the fixed all-wrong comparator, so a calibrated cell-matching attribution is supported. If subsequent checks were to make that contrast absent or mixed, attribution to cell matching should be removed.

This guide identifies evidence-supported framing only; it does not implement a manuscript revision.
