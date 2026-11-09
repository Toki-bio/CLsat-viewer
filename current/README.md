# Current Darevskia catalogue viewer (revision, 2026-11-09)

Two viewers, side by side:

- **Old** (as submitted, before the review): https://toki-bio.github.io/CLsat-viewer/index.html
- **Current** (this folder): https://toki-bio.github.io/CLsat-viewer/current/index.html

The current viewer is built by `satpipe viewer` (Toki-bio/CLsat_workflow main b08df72) from run
`monsoon:/scratch/sk4386/satpipe_final/run_z3_v2`: the final catalogue (run_z3) with the
haplotype-duplicate rule and the stage-19 orthology rules.

- 550 loci (run_z3 had 555). Removed as haplotype duplicates of *unim*: unif L34, L45, L46, L48, L49.
- 19 ortholog locus pairs (stage 19: tandem and repeat masking, single-copy flanks, clear best 0.9,
  reciprocal, flanks extended 600/3000/5000 when unresolved).
- 4 candidate loci outside the catalogue.

Loci whose family call changed against the old viewer (window majority agrees with the new call):
L389.2 (199/200 windows CLsat1), L51.2 (194/200 CLsat1), L4 (193/200 CLsat1), L318.2 (199/200 CLsat3).
