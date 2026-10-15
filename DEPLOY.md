# Deploying this viewer — read before changing anything

There is **one** live page and **one** deploy target. Everything else in this repo, and
every other clone on disk, is history.

## The live page

<https://toki-bio.github.io/CLsat-viewer/> — served from the **root** of this repository.
There is no other page.

`corrected/` no longer holds a viewer. From 2026-07-29 it contains a single redirect to
the root, so that any previously shared `/corrected/` link still resolves. Do not put data
files back into it, and do not deploy there. The split between a stale root page and a
maintained `corrected/` subpage was the single worst source of confusion in this project's
history — the root page served a superseded 706-locus catalogue for months while corrected
work accumulated out of sight.

## Deploy target

This clone: `C:\Users\T\Documents\clsat_viewer_update\repo_deploy2` (remote
`github.com/Toki-bio/CLsat-viewer`, branch `main`).

Four other clones exist on disk and are all stale — two from February 2026 at v2.10, two
from June 2026 whose remote URL embeds a credential (revoke that token on GitHub and delete
those clones). Do not use them. If you find yourself in one, stop.

## Where the current files are built

`C:\work\CLsat\viewer_v2\` is the working directory. Deploy = copy its runtime set to this
repo root, commit, push.

Runtime set (everything else in `viewer_v2` is helper scripts and backups — do not copy):

```
index.html
flanks.json  subfamily.json  flank_orthologs.json
full_flank_ortho_table.tsv  scaffold_lengths.tsv
MASTER.merged.bed  MASTER.final_assignment.clean.tsv  MASTER.locus_bits.tsv
MASTER_flank.bits.tsv.gz
MASTER_CLsat1.bits.tsv.gz  MASTER_CLsat2.bits.tsv.gz  MASTER_CLsat3.bits.tsv.gz  MASTER_CLsat4.bits.tsv.gz
CLsat1.cons.fa  CLsat2.cons.fa  CLsat3.cons.fa  CLsat4.cons.fa
```

Since v2.25 the five per-window bitscore tracks are committed **gzip-compressed** (`gzip -9 -n`, 11 MB
instead of 161 MB; the page inflates them with `DecompressionStream`, and falls back to a plain
`.tsv` of the same name if the `.gz` is missing). Do not commit the uncompressed files.

`terminal_monomers.json` is optional. The file that shipped until 2026-10-15 was built from the
superseded 707-locus catalogue: 260 current loci were missing and 115 shared IDs carried another
locus's record. It is kept as `terminal_monomers.707locus.withdrawn.json` and is **not loaded**.
Regenerate it with the corrected `05_subfamily_terminal_boundary/terminal_monomers.sh` from the
workflow repo, then save it as `terminal_monomers.json`; the page shows the panel only for records
whose monomer length matches the locus family. `boundary_phase*.json` and `terminal_hist.json`
were empty or unused and are gone.

Previous pages live in git history (`git log -- index.html`); do not keep `index_*.html` archives
in the Pages root, because they are served live.

Leave `README.md` and `VIEWER_MANUAL.md` alone unless the code they describe changed.

## Δ sign convention in `MASTER.final_assignment.clean.tsv`

Column 6 is **assigned-family bits − strongest other family**, for every flag (since 2026-10-15).
Before that date the `UNIQUE_SEED_*` rows carried the opposite sign (best other − seed), so the
viewer printed a mixed column. `arbitration_patched.sh` in the workflow repo now writes the uniform
convention; a TSV produced by an older run must be sign-flipped on those rows before deploy.

## Version

`index.html` contains `VIEWER_VERSION`. **Bump it on every deploy.** If the local and live
copies report the same version while differing in content, the string is worthless — that
happened between v2.19 and this deploy.

## History of what the root page served

| version | catalogue | note |
|---|---|---|
| unversioned (to 2026-07-29) | **706 loci** | superseded; pre-correction data. Archived as `index_archive_706locus_pre2026-07-29.html` |
| v2.19-clip-20260727 | 627 loci | only ever at `corrected/`, never at root |
| v2.20-candidates-20260728 | 634 loci | 627 catalogue + 7 uncatalogued candidates (flag `CANDIDATE_UNCATALOGUED`); `corrected/` only |
| **v2.21-overlaid-20260729** | 634 loci | **first version served at the root.** Adds the "bars are overlaid, not stacked" banner and the 1 kb segment-boundary tolerance note; `corrected/` reduced to a redirect |
| v2.24-loader-20260729 | 634 loci | browser caching allowed; `getBaseUrl` hardened |
| v2.25-audit-20261015 | 634 loci | audit fixes: Smith–Waterman in a Web Worker with 3 bytes/cell and single-pass hit extraction; gzip tracks; L541 in the full table; uniform Δ sign; complete flag vocabulary; stale terminal-monomer file withdrawn; manuals rewritten from current data |

The root page served the superseded 706-locus catalogue for months while corrected work
accumulated in `corrected/`. That is the mistake this file exists to prevent.

## Catalogue counts, for checking a deploy landed

`wc -l MASTER.merged.bed` at the root should equal the locus count of the version you
just shipped (634 for v2.20 onward: 627 catalogue loci + C1–C7). If it says 706, the deploy did not happen.
The embedded lists must match too: `FULL_TABLE` 627 rows, `FLAGGED_LOCI` 35 (28 + 7 candidates),
`ORTHO_LOCI` 91.
