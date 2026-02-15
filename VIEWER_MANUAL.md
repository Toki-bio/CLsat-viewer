# CLsat Locus Detail Viewer – Manual

## Overview

`locus_detail_viewer.html` is a two-panel interactive viewer for inspecting
CLsat satellite repeat loci across multiple species and scaffolds.

- **Top panel (Scaffold Overview):** shows all loci on a scaffold with overlaid
  per-window bitscore bars (CLsat1–4), colored by consensus. Click a locus
  rectangle to open it in the detail panel.
- **Bottom panel (Per-window Detail):** shows individual scan-window bitscores
  for the selected locus, with zoom/pan. Supports BED coordinate export,
  NCBI FASTA fetching, and monomer analysis via a local server.

---

## How to launch

### Inside VS Code (Five Server)
The workspace root must be served. Navigate to:
```
http://127.0.0.1:5500/viewer/files/out/locus_detail_viewer.html
```
(Five Server serves from the workspace root on port 5500 by default.)

### Outside VS Code (standalone)
Double-click **`open_viewer.bat`** in `viewer/files/out/`.
It starts a PowerShell HTTP server on port 8080 and opens the viewer
in your default browser. Keep the console window open; close it to stop.

No Python or Node.js required — uses PowerShell's built-in .NET HttpListener.

---

## Files in `viewer/files/out/`

### The viewer
| File | Description |
|------|-------------|
| `locus_detail_viewer.html` | Main viewer (current version, ~1360 lines) |
| `open_viewer.bat` | One-click launcher for standalone use |
| `serve.ps1` | PowerShell HTTP server script used by the launcher |

### Historical viewer versions (read-only reference)
| File | Description |
|------|-------------|
| `locus_detail_viewer_stable.html` | Early stable snapshot |
| `locus_detail_viewer_before_fasta.html` | Before NCBI FASTA integration |
| `locus_detail_viewer_with_ncbi_fasta.html` | With FASTA, before monomer server |
| `locus_viewer.html` | Predecessor (no per-window detail panel) |
| `test.html` | Earlier MASTER viewer (overview only, no detail) |
| `clsat_bits_viewer.html` | Simple single-file bitscore viewer |

### Required data files (auto-loaded on startup)
| File | Format | Description |
|------|--------|-------------|
| `MASTER.merged.bed` | BED (scf, start, end, name) | Locus coordinates. Locus IDs (L1, L2, …) generated from line number. The `name` column encodes `species\|scf:start-end\|MASTER`. |
| `MASTER.final_assignment.clean.tsv` | TSV with header | Seed-priority arbitration results. Columns: `locus_id, seed_subfamily, seed_bits, best_other, best_other_bits, delta_bits, FLAG`. FLAG is one of: `SEED_ONLY`, `SEED_WITH_COMPETITOR`, `SEED_CONTRADICTION`. |
| `MASTER.locus_bits.tsv` | TSV (no header) | Per-locus max bitscore summary. Columns: `locus_id, CLsat1_max, CLsat2_max, CLsat3_max, CLsat4_max`. |
| `MASTER_CLsat1.bits.tsv` | TSV with header | Per-window scan bitscores for CLsat1. Columns: `scaffold, pos, bits, name, species`. Window step = **145 bp**. |
| `MASTER_CLsat2.bits.tsv` | TSV with header | Same for CLsat2. Window step = **145 bp**. |
| `MASTER_CLsat3.bits.tsv` | TSV with header | Same for CLsat3. Window step = **146 bp**. |
| `MASTER_CLsat4.bits.tsv` | TSV with header | Same for CLsat4. Window step = **116 bp**. |

**Important:** Each CLsat has a different monomer length and therefore a different
scan step. The steps are determined by the consensus length minus 1:

| Consensus | Length | Scan step |
|-----------|--------|-----------|
| CLsat1 | 146 bp | 145 bp |
| CLsat2 | 145 bp | 145 bp |
| CLsat3 | 146 bp | 146 bp |
| CLsat4 | 116 bp | 116 bp |

The consensus FASTA files are also present in the folder:
`CLsat1.cons.fa`, `CLsat2.cons.fa`, `CLsat3.cons.fa`, `CLsat4.cons.fa`.

### Optional: Monomer analysis server
| File | Description |
|------|-------------|
| `monomer_server.py` | FastAPI server for monomer decomposition (port 8000) |
| `server_control.py` | Control server to start/stop monomer_server (port 9001) |
| `start_monomer_server.bat` | Launcher for monomer_server.py |
| `start_server_control.bat` | Launcher for server_control.py |
| `start_viewer.bat` | Combined launcher (servers + viewer) |

**Note:** These require a working Python 3.10+ with FastAPI/uvicorn.
The `.venv/` in this folder has a broken symlink to a deleted Python 3.10
install. To fix, delete `.venv/` and recreate:
```
python -m venv .venv
.venv\Scripts\activate
pip install fastapi uvicorn
```

---

## Interaction guide

### Navigation
- **Species dropdown** → **Scaffold dropdown** → overview draws
- **Jump dropdown** lists contradictions (⚠) and competitors (●)
- **Find locus** text input: type e.g. `L144` and click Go
- **Reset zoom** resets both panels

### Overview panel
- **Scroll-zoom** to zoom, **drag** to pan
- **Click a locus rectangle** to select it and show per-window detail
- Locus outlines: white border = selected, thin = unselected
- Bar colors: blue = CLsat1, green = CLsat2, pink = CLsat3, yellow = CLsat4

### Detail panel
- **Scroll-zoom** and **drag** to pan
- **BED coordinates** update as you pan/zoom (shows visible region)
- **Copy BED** copies `scaffold\tstart\tend` to clipboard
- **Get FASTA** fetches the visible region sequence from NCBI E-utilities
- **Analyze Monomers** sends the FASTA to monomer_server.py for decomposition

### Filters
- `ALL` / `SEED_ONLY` / `SEED_WITH_COMPETITOR` / `SEED_CONTRADICTION`
  filters loci shown in the overview

---

## Known issues and limitations

### 1. Monomer analysis does NOT replicate searsat16

The "Analyze Monomers" button calls `monomer_server.py` which uses a
**sliding-window identity scan** (`find_monomers()`) and a custom
**Python DP alignment** (`simple_align()`).

The original `searsat16_reference.sh` pipeline does:
1. **ssearch36** Smith-Waterman alignment (`ssearch36 -r +2/-2 -g -3 …`)
   to find monomer matches — a proper gapped local alignment, not a
   fixed-length window identity check
2. **bedtools merge** with strand-awareness to merge overlapping hits
3. **Gap-based array splitting** (>100bp gaps = separate arrays)
4. **bedtools getfasta** to extract array sequences with ±100bp flanks
5. **MAFFT** multiple sequence alignment (`--localpair --maxiterate 1000`)
   of consensus + extracted array pieces

The monomer_server.py replacement:
- Uses a **fixed-length sliding window** instead of ssearch36's gapped
  Smith-Waterman — this will miss monomers with insertions/deletions
  and truncated terminal copies
- Uses a **simple O(mn) Python DP pairwise alignment** instead of MAFFT
  MSA — this is extremely slow for long loci and produces a 2-sequence
  alignment rather than a proper multi-monomer MSA
- Only has **CLsat1 consensus hardcoded** in `DEFAULT_CONSENSUSES` —
  CLsat2/3/4 are missing, so analysis fails for non-CLsat1 loci unless
  `consensus_seq` is passed
- The alignment output format (repeated consensus vs locus) doesn't
  match searsat16's individual-monomer MAFFT output

**Status:** This feature is non-functional without ssearch36 and MAFFT
being available on the serving machine. The Python DP approach is not
a viable substitute.

### 2. NCBI FASTA fetch requires valid accessions
The "Get FASTA" button calls NCBI E-utilities with the scaffold name as
an accession. This only works if your scaffolds are deposited NCBI
accessions (e.g., JAWWNF010000123.1). It will fail for custom/local
scaffold names.
