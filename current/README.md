# Current Darevskia catalogue viewer (revision, 2026-11-10)

Two viewers, side by side:

- **Old** (as submitted, before the review): https://toki-bio.github.io/CLsat-viewer/index.html
- **Current** (this folder): https://toki-bio.github.io/CLsat-viewer/current/index.html

The current viewer is built by `satpipe viewer` (Toki-bio/CLsat_workflow main 842945a) from run
`monsoon:/scratch/sk4386/satpipe_final/run_gb_v2`, the final catalogue of the paper:

- Six genome copies: *D. valentini* (val), *D. mixta* (mix), the two parent-split copies of *D. armeniaca*
  (armf, *D. valentini* side; armm, *D. mixta* side) and the GenBank *D. unisexualis* pair, labelled by parent:
  unv = GCA_032766585.1 (GenBank "mHiFi", *D. valentini* side), unn = GCA_032766555.1 (GenBank "pHiFi",
  *D. nairensis* side).
- 594 loci. The haplotype-duplicate rule removed unv L134 (47,681 bp, CLsat1.2), an exact duplicate of part of
  unn L67. unv L133 (84 kb, CLsat1.2) is identical to unn L80, but both sit on flankless contigs, so it is kept
  and flagged (probably one array assembled in both copies).
- 19 ortholog locus pairs from the stage-19 flank rules (tandem and repeat masking, single-copy flanks, clear best
  0.9, reciprocal, flanks extended 600/3000/5000 when unresolved). Orthology in the paper rests on HiFi reads.
- 4 candidate loci outside the catalogue.

The previous build of this folder (run_z3_v2, trio *D. unisexualis* copies unif/unim, 550 loci) is in the git
history (commit 208fe0e).
