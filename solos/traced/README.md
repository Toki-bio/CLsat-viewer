# The 29 dispersed CLsat solo copies traced to the other five assemblies

`trace.sh` (run on KIT, 2026-11-18, genomes = the six catalog assemblies, checked by sequence). Per copy: on each side
the nearest 400-bp window (0-4.6 kb from the copy) that is single-copy in its own assembly (one megablast hit >= 200 bp
at >= 90%); each window searched in the other five assemblies (best hit, second best < 0.9 x best bits, hit must reach the
window end next to the copy); a place = both windows on one contig, same orientation, in order, <= 50 kb apart.
Between them, CLsat by ssearch36 against the doubled family consensus (>= 60 bits).

`status.tsv`: id, source assembly, sub-subfamily, copy length, target, left hit, right hit, status (present, empty,
other_sequence, one_flank_only, not_found, flanks_disagree), gap between the projected copy ends (bp), CLsat bp in the gap.
`<id>.aln.fa`: MAFFT L-INS-i of the source copy with 300 bp of flank and every target place (both flanks; or one flank
plus copy length plus 300 bp, marked one_flank_only). Controls: positions P01 and P07 of the catalog's stage-19
orthology are recovered.
