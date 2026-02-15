# CLsat Locus Detail Viewer

Interactive two-panel viewer for inspecting centromeric satellite repeat (CLsat) loci
across multiple species and scaffolds. Built with D3.js — runs entirely in the browser
with no server dependencies.

## Features

- **Scaffold Overview** — all loci on a scaffold with per-window bitscore bars (CLsat1–4)
- **Per-window Detail** — zoom into individual loci, see bitscore patterns at each scan window
- **NCBI FASTA fetch** — retrieve locus sequences directly from NCBI E-utilities
- **Client-side monomer analysis** — Smith-Waterman alignment replicating the
  [searsat16](searsat16_reference.sh) pipeline entirely in JavaScript (no server needed)
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
> over HTTP — opening the HTML file directly (`file://`) won't work.

## Data Files

| File | Description |
|------|-------------|
| `MASTER.merged.bed` | Locus coordinates (BED format) |
| `MASTER.final_assignment.clean.tsv` | Seed-priority arbitration results |
| `MASTER.locus_bits.tsv` | Per-locus max bitscore summary |
| `MASTER_CLsat{1-4}.bits.tsv` | Per-window scan bitscores for each consensus |
| `CLsat{1-4}.cons.fa` | Consensus sequences |

## Monomer Analysis

The "Analyze Monomers" button runs a client-side Smith-Waterman local alignment
of the selected CLsat consensus against the fetched locus sequence. The algorithm
replicates the key steps of `searsat16_reference.sh`:

1. **Smith-Waterman** with ssearch36 scoring (match +2, mismatch −2, gap open −16, gap extend −3)
2. **Merge** overlapping hits (strand-aware, like `bedtools merge -s -d -1`)
3. **Filter** by identity (≥75%) and length (≥90% of consensus)
4. **Split** into tandem arrays at gaps >100 bp
5. **Require** ≥3 monomers per array

## Documentation

See [VIEWER_MANUAL.md](VIEWER_MANUAL.md) for detailed file descriptions, interaction
guide, and known limitations.

## License

This project is provided for academic and research use.
