# CLsat, Agi160 and the 147-bp arrays of Lacerta agilis

- `units_CLsat_Agi160_147bp.aln.fa`: the nine CLsat sub-subfamily consensuses, the unit of the two 147-bp arrays of
  L. agilis (NC_046330.1:466006-467827 and 479892-481857, GCF_009819535.1; both arrays have the same TRF unit), the 18
  Agi160 GenBank clones (AY184824-35, AY652647-52) and the TRF units (period 140-200 bp) of the L. agilis Agi160 arrays.
- `Lagilis_147bp_arrays_monomers.aln.fa`: CLsat3.2, the 147-bp unit and the monomers of both arrays, cut in the unit's frame.

Common frame before alignment (`prep_common_frame.py`): CLsat consensuses as deposited; the 147-bp unit rotated and
stranded onto CLsat3.2; every Agi160 sequence onto clone AY184829, itself framed on the 147-bp unit. Alignment:
MAFFT 7 L-INS-i (`--localpair --maxiterate 1000`). TRF records from the outgroup scan of the CLsat revision
(TRF 2 5 7 80 10 50 2000).

2026-11-18: flanks added to `Lagilis_147bp_arrays_monomers.aln.fa`. The region NC_046330.1:464,006-483,857 (chromosome W
of L. agilis, rLacAgi1.pri) was fetched from NCBI (`NC_046330_464006_483857.fa`); both TRF arrays match it base for base.
Rows `<array>_left_flank_N` / `_right_flank_N` (lowercase) hold the genome sequence next to each array, continued in the
frame of the 147-bp unit without alignment, so the base next to the array sits in the column next to the array's end.
