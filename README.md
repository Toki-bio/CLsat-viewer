# CLsat Locus Detail Viewer

Interactive two-panel viewer for inspecting centromeric satellite repeat (CLsat) loci
across multiple species and scaffolds. Built with D3.js — runs entirely in the browser
with no server dependencies.

## Features

- **Scaffold Overview** — all loci on a scaffold with per-window bitscore bars (CLsat1–4)
- **Per-window Detail** — zoom into individual loci, see bitscore patterns at each scan window
- **NCBI FASTA fetch** — retrieve locus sequences directly from NCBI E-utilities
- **Client-side monomer analysis** — Smith-Waterman alignment approximating the
  [searsat16](searsat16_reference.sh) pipeline entirely in JavaScript, in a Web Worker (no server needed)
- Locus filtering by arbitration flag (SEED_ONLY / SEED_WITH_COMPETITOR / SEED_CONTRADICTION)
- BED coordinate export, zoom/pan, jump to contradictions

## Quick Start

Open `index.html` in any modern browser via a local HTTP server:

```bash
# Python
python -m http.server 8080

# Node.js
npx serve .

# Or use the GitHub Pages deployment directly
```

> **Note:** The viewer uses `fetch()` to load data files, so it must be served
> over HTTP — opening the HTML file directly (`file://`) won't work. Locus sequences
> are fetched live from NCBI E-utilities, so an internet connection is needed for
> **Get FASTA** and **Analyze Monomers**.

## Data Files

| File | Description |
|------|-------------|
| `MASTER.merged.bed` | Locus coordinates (BED, 0-based half-open); 627 catalogue loci + 7 uncatalogued candidates C1–C7 |
| `MASTER.final_assignment.clean.tsv` | Arbitration results, no header: id, family, bits, strongest other family, its bits, Δ (= col 3 − col 5), flag, seed families |
| `MASTER.locus_bits.tsv` | Per-locus max bitscore per family |
| `MASTER_CLsat{1-4}.bits.tsv.gz` | Per-window scan bitscores for each consensus (gzip; inflated in the browser) |
| `MASTER_flank.bits.tsv.gz` | Per-window bitscores in the flanks, for the noise floor |
| `flanks.json`, `flank_orthologs.json`, `subfamily.json` | Flank classes, cross-assembly orthologs, sub-subfamily calls |
| `scaffold_lengths.tsv` | Scaffold lengths, for contig-edge markers (5 of 406 scaffolds absent) |
| `full_flank_ortho_table.tsv` | Source table for the embedded flank-orthology list |
| `CLsat{1-4}.cons.fa` | Consensus sequences (also embedded in `index.html`) |
| `terminal_monomers.json` | Optional; see DEPLOY.md — the previous file was withdrawn as stale |

## Monomer Analysis

The "Analyze Monomers" button runs a client-side Smith-Waterman local alignment
of all four CLsat consensuses against the fetched locus sequence (both strands), in a
Web Worker. It follows `searsat16_reference.sh` in these steps:

1. **Smith-Waterman** with ssearch36 scoring (match +2, mismatch −2, first gap residue −16, gap extend −3)
2. **Merge** overlapping same-strand hits (like `bedtools merge -s -d -1`), then
3. **Filter** by identity (≥75 %, adjustable) and length (≥90 % of consensus)
4. **Split** into tandem arrays at true gaps >100 bp

It is an approximation, not a re-implementation. Known differences from the shell pipeline:

- Hits are accepted above a fixed score (35 % of a perfect monomer match, ≈102 bits for CLsat1)
  instead of an E-value cut-off, and `N` in a consensus scores as a mismatch.
- Minus-strand hits are kept and shown as minus-strand arrays; `searsat16` discards them.
- All four consensuses are scanned; overlapping hits are resolved by score and each monomer is
  then re-scored against every consensus, so a locus can show a CLsat1→CLsat3 transition.
- Arrays need ≥3 monomers (`MIN_MONOMERS_PER_ARRAY`); `searsat16` needs 4. After a family
  transition, interior sub-arrays need ≥2 and edge sub-arrays ≥1.
- Merged hits are capped at 1.3× the consensus length and wider hits are decomposed.
- No ±100 bp array expansion and no MAFFT: the "MSA" is a consensus-anchored stack of pairwise
  alignments (insertions relative to the consensus are not displayed). Use **Copy All FASTA** to
  align the monomers with MAFFT yourself.
- The browser scan is limited to 200 kb of sequence (`SW_MAX_BP`).

## Documentation

See [VIEWER_MANUAL.md](VIEWER_MANUAL.md) for detailed file descriptions, interaction
guide, and known limitations.

## License

This project is provided for academic and research use.
