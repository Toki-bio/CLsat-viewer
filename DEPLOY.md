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
from June 2026 with a personal access token embedded in the remote URL. Do not use them.
If you find yourself in one, stop.

## Where the current files are built

`C:\work\CLsat\viewer_v2\` is the working directory. Deploy = copy its runtime set to this
repo root, commit, push.

Runtime set (everything else in `viewer_v2` is helper scripts and backups — do not copy):

```
index.html
flanks.json  subfamily.json  flank_orthologs.json
boundary_phase.json  boundary_phase_hist.json
full_flank_ortho_table.tsv  scaffold_lengths.tsv
MASTER.merged.bed  MASTER.final_assignment.clean.tsv  MASTER.locus_bits.tsv
MASTER_flank.bits.tsv
MASTER_CLsat1.bits.tsv  MASTER_CLsat2.bits.tsv  MASTER_CLsat3.bits.tsv  MASTER_CLsat4.bits.tsv
CLsat1.cons.fa  CLsat2.cons.fa  CLsat3.cons.fa  CLsat4.cons.fa
```

Before overwriting, archive the outgoing page as `index_pre_<version>_<date>.html`.

Leave `README.md`, `VIEWER_MANUAL.md` and the `index_*.html` archives alone.

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

The root page served the superseded 706-locus catalogue for months while corrected work
accumulated in `corrected/`. That is the mistake this file exists to prevent.

## Catalogue counts, for checking a deploy landed

`wc -l MASTER.merged.bed` at the root should equal the locus count of the version you
just shipped (634 for v2.20). If it says 706, the deploy did not happen.
