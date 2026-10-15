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
changes along the array**. The catalogue holds **627** such loci (540 MASTER loci before segmentation; 626 from the pipeline plus L541, added by hand after a singleton screen). Seven further CLsat-matching arrays found by that screen but deliberately left out of the catalogue are shown as `C1`–`C7` for inspection only. IDs with a suffix (`L36.1`, `L36.2`) are segments
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

- **Monomer list** — every monomer with coordinates, length and % identity to the consensus (from the
  Smith–Waterman traceback, not from the scan bitscore).
- **Monomer map** — the array drawn as arrows; direction = monomer strand. Hover to highlight.
- **MSA** — a consensus-anchored stack: each monomer is aligned pairwise to the consensus and projected onto
  consensus columns (insertions relative to the consensus are not shown), with up to 300 bp of genomic flank attached
  to the first and last row. It is not a MAFFT multiple alignment; use **Copy All FASTA** to run one.

## Terminal monomers vs reference

Displayed **below** the monomer analysis. It answers one narrow question: does the array *begin* and
*end* on a whole monomer, or is the first/last monomer only partial?

- Choose the reference: **sub-subfamily**, **subfamily**, or **this-locus** consensus.
- `ref` row — the reference consensus. Dim grey = positions this monomer does **not** cover
  (missing coverage, *not* extra bases).
- `mono` row — the actual genomic monomer. Underlined = matches the reference; red background = mismatch.
- `trim5` / `trim3` — how many reference bases are missing from the start / end. These are descriptive: both zero is
  consistent with a complete terminal monomer under the current rotation convention, but the panel makes no complete/partial call.
- Dim lowercase at the outer edge — 5 bp of real genomic flank *outside* the array. If the array runs
  into a contig edge there is no flank, and the header says *(contig edge)*.

> **Caveat.** A tandem monomer has no intrinsic start. "Complete" and "partial" are defined relative to
> a rotation convention inherited from the original published CLsat clones — not a biological landmark.
> A different, equally defensible phase would relabel some termini without changing any sequence.

---

## Inspection lists

Three curated sets, opened from the sidebar buttons (orthologous, problematic, and the full flank-orthology table). In all of them, **clicking any row opens that locus**.

### Orthologous loci

All **91** loci with at least one usable unique flank, scored by their positional orthologs in other assemblies
(flank identity ≥80 %, ≥200 bp, **direct partners only** — transitive chains are not followed). A locus can have
several partner loci in the same assembly; in this dataset no locus has partners in more than two other assemblies.

| Preset | Loci |
|---|---|
| Has an ortholog | 58 |
| In both other lineages (partners in 2 assemblies) | 16 |
| Private (no cross-assembly ortholog) | 33 |
| Strictly private (both flanks unique, still no ortholog) | 2 |

The partner chips are clickable: use them to **jump straight to the same genomic position in another
genome** and see what satellite sits there.

Two things this list makes visible:

- **Conserved position ≠ conserved content.** 23 of the 53 orthologous locus pairs (43 %) carry a *different* CLsat
  family at the two ends.
- **Anchorability is uneven.** Of the 91 anchorable loci, 30 are in *D. mixta*, 22 in `unp`, 19 in *D. nairensis*;
  the 33 private loci are 18 CLsat1, 14 CLsat2 and 1 CLsat3.

### Problematic loci

The **28** catalogue loci flagged by arbitration (24 segment calls and 2 seed-conflict calls with a margin under 15 bits,
L186, and L540 which lacked a seed), plus the 7 uncatalogued candidates. Columns show the assigned call and bitscore,
the strongest other family and bitscore, and the margin Δ = assigned − strongest other (negative means the other family
scored higher).

Start with the single `UNIQUE_SEED_CONTRADICTION` locus (**L186**) — the only case where the independent rescan
confidently disagrees with the seed (CLsat1 by 31 bits over the seed call CLsat3, which is retained).

---

## Arbitration flags

| Flag | Meaning |
|---|---|
| `UNIQUE_SEED_ONLY` | One seed family; no other family scored in the rescan. |
| `UNIQUE_SEED_WITH_COMPETITOR` | Clean seed call; a competitor scores but is within 20 bits. |
| `UNIQUE_SEED_CONTRADICTION` | The rescan favours a **different** family by ≥ 20 bits; the seed call is retained and the locus flagged. Only L186. |
| `SEED_CONFLICT_RESOLVED` | Seeds from more than one subfamily fell in the locus; resolved by bitscore. |
| `SEED_CONFLICT_RESOLVED_AMBIGUOUS` | As above, but the margin is under 15 bits. |
| `SEGMENT_RESOLVED` | Locus split at a sustained subfamily transition; each segment called cleanly. |
| `SEGMENT_RESOLVED_AMBIGUOUS` | Split locus whose segment margin is under 15 bits. |
| `NO_SEED_RESOLVED` | No seed fell in the locus; family from the rescan alone (L540). |
| `CANDIDATE_UNCATALOGUED` | Not a catalogue member; one of C1–C7. |

---

## Data files

The viewer is a single self-contained `index.html` that fetches static products sitting beside it:

`MASTER.merged.bed` · `MASTER.final_assignment.clean.tsv` · `MASTER.locus_bits.tsv` ·
`MASTER_CLsat{1..4}.bits.tsv.gz` · `MASTER_flank.bits.tsv.gz` · `flanks.json` · `flank_orthologs.json` · `subfamily.json` ·
`scaffold_lengths.tsv` · `CLsat{1..4}.cons.fa`

The per-window bitscore files are gzip-compressed; the page decompresses them in the browser.
`terminal_monomers.json` is loaded when present; the shipped file was withdrawn on 2026-10-15 because it predated the
catalogue correction (see DEPLOY.md), so the terminal-monomer panel is hidden until it is regenerated.

The three inspection lists are **embedded in `index.html`** and need no extra files.

## Tips

- **Esc** closes any panel; every table has a search box.
- If **Analyze Monomers** looks wrong, check that you fetched the FASTA first and that the locus is zoomed in.
- The per-window bitscore track is the ground truth — a table label can hide a locus whose signal switches half-way.
