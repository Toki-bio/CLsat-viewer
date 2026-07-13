# CLsat Locus Viewer — Manual

An interactive browser for the CLsat satellite arrays catalogued across six *Darevskia* genome
assemblies. It exists so that **any locus-level claim can be checked against the raw evidence** —
per-window bitscores, monomer structure, array boundaries and cross-assembly orthology — rather than
being taken on trust from a summary table.

The same text is available inside the app: click **? Manual** next to *Scaffold Overview*.

---

## Assemblies

| Code | Species |
|------|---------|
| `mix` | *D. mixta* |
| `arm` | *D. armeniaca* |
| `val` | *D. valentini* |
| `unp` | *D. unisexualis*, paternal-haplotype assembly |
| `unm` | *D. unisexualis*, maternal-haplotype assembly |
| `nai` | *D. nairensis* |

## What a "locus" is

A **segmented MASTER locus**: a contiguous CLsat-similar domain, expanded outward from a TRF seed to
its true array boundary, merged with overlapping domains, and then **split where the dominant subfamily
changes along the array**. There are 707 such loci. IDs with a suffix (`L36.1`, `L36.2`) are segments
of one longer array that contained a sustained subfamily transition.

---

## Quick workflow

1. **Select a locus** — click a row in the **Locus Table**, click an array in the **Scaffold Overview**,
   or type an ID (e.g. `L152`) into **Find locus** and press **Go**.
2. **Zoom in** — mouse-wheel to zoom, drag to pan. Zoom until individual monomers resolve; the
   monomer-level tools need the array readable, not the whole scaffold.
3. **Get FASTA** — fetches the locus sequence.
4. **Analyze Monomers** — splits that sequence into individual monomers.
5. **Inspect below** — monomer list → monomer map → MSA → terminal-monomer trim.

---

## Panels

- **Scaffold Overview** — every locus on the selected scaffold; bar height = bitscore.
- **Per-window Bitscores** — the selected locus scanned against all four CLsat consensuses, one bar per
  scan window. This is *why* a locus got its call. A composite array shows the dominant colour switching
  part-way along; a contested locus shows two colours running neck-and-neck. **When a table label and
  this track disagree, believe this track.**
- **Locus Table** — filter by ID, assembly, family or flag. The *Show* and *Subtype* filters stay in
  sync with the overview.
- **Locus Ends & Flanks** — the unique sequence flanking the array, its truncation status, and the
  cross-assembly flank matches that define positional orthology.

## Monomer analysis

- **Monomer list** — every monomer with coordinates, length and % identity to the consensus.
- **Monomer map** — the array drawn as arrows; direction = monomer strand. Hover to highlight.
- **MSA** — monomers stacked and aligned with their genomic flanks attached, so substitutions that
  recur down the array (candidate higher-order structure) appear as vertical stripes.

## Terminal monomers vs reference

Displayed **below** the monomer analysis. It answers one narrow question: does the array *begin* and
*end* on a whole monomer, or is the first/last monomer only partial?

- Choose the reference: **sub-subfamily**, **subfamily**, or **this-locus** consensus.
- `ref` row — the reference consensus. Dim grey = positions this monomer does **not** cover
  (missing coverage, *not* extra bases).
- `mono` row — the actual genomic monomer. Underlined = matches the reference; red background = mismatch.
- `trim5` / `trim3` — how many reference bases are missing from the start / end. Both zero means the
  terminal monomer is complete under the current rotation convention.
- Dim lowercase at the outer edge — 5 bp of real genomic flank *outside* the array. If the array runs
  into a contig edge there is no flank, and the header says *(contig edge)*.

> **Caveat.** A tandem monomer has no intrinsic start. "Complete" and "partial" are defined relative to
> a rotation convention inherited from the original published CLsat clones — not a biological landmark.
> A different, equally defensible phase would relabel some termini without changing any sequence.

---

## Inspection lists

Two curated sets, opened from the sidebar buttons. In both, **clicking any row opens that locus**.

### Orthologous loci

All **550** loci with at least one usable unique flank, scored by how many *other* assemblies they are
positionally orthologous to (flank identity ≥80 %, ≥200 bp, **direct partners only** — transitive chains
are not followed, because the orthology graph percolates into one giant component).

| Preset | Loci |
|---|---|
| Conserved in all 5 other assemblies | 74 |
| In ≥3 assemblies | 296 |
| Private (no cross-assembly ortholog) | 93 |
| Strictly private (both flanks unique, still no ortholog) | 22 |

The partner chips are clickable: use them to **jump straight to the same genomic position in another
genome** and see what satellite sits there.

Two things this list makes visible:

- **Conserved position ≠ conserved content.** 42 % of orthologous locus pairs (614 / 1,461) carry a
  *different* CLsat subfamily at the two ends — most often CLsat1 opposite CLsat3.
- **Positional novelty is lineage-specific.** The strictly private loci concentrate in *D. mixta*
  (14 of its 56 both-flank-anchorable loci) and in CLsat2, while *D. nairensis* — despite having the most
  CLsat loci of any assembly — has none.

### Problematic loci

The **27** loci flagged by arbitration: composite arrays, short residual segments, and near-tied family
calls. Columns show the seed call and bitscore, the strongest competitor and bitscore, and the margin Δ.

Start with the two `UNIQUE_SEED_CONTRADICTION` loci (**L191**, **L425**) — those are the only cases where
the independent rescan confidently disagrees with the seed.

---

## Arbitration flags

| Flag | Meaning |
|---|---|
| `UNIQUE_SEED_WITH_COMPETITOR` | Clean seed call; a competitor scores but does not threaten it. |
| `UNIQUE_SEED_CONTRADICTION` | The rescan confidently favours a **different** family than the seed (Δ ≥ 20 bits). Only 2 loci. |
| `SEED_CONFLICT_RESOLVED` | Seeds from more than one subfamily fell in the locus; resolved by bitscore. |
| `SEED_CONFLICT_RESOLVED_AMBIGUOUS` | As above, but the margin is small. |
| `SEGMENT_RESOLVED` | Locus split at a sustained subfamily transition; each segment called cleanly. |
| `SEGMENT_RESOLVED_AMBIGUOUS` | Split locus whose segment boundary or call is uncertain. |

---

## Data files

The viewer is a single self-contained `index.html` that fetches static products sitting beside it:

`MASTER.merged.bed` · `MASTER.final_assignment.clean.tsv` · `MASTER.locus_bits.tsv` ·
`MASTER_CLsat{1..4}.bits.tsv` · `flanks.json` · `flank_orthologs.json` · `subfamily.json` ·
`terminal_monomers.json` · `terminal_hist.json` · `boundary_phase.json` · `boundary_phase_hist.json` ·
`CLsat{1..4}.cons.fa`

The two inspection lists are **embedded in `index.html`** and need no extra files.

## Tips

- **Esc** closes any panel; every table has a search box.
- If **Analyze Monomers** looks wrong, check that you fetched the FASTA first and that the locus is zoomed in.
- The per-window bitscore track is the ground truth — a table label can hide a locus whose signal switches half-way.
