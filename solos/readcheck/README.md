# HiFi read check of the traced solo places (2026-11-18, KIT)

`readcheck.sh`: for the 12 dispersed copies with single-copy flank windows on both sides and at least one place found
in another assembly (`ids_tested.txt`), reads of each assembled animal holding window words were fished from the HiFi BAMs
(D. valentini 218 Dv3; D. mixta 80 Dm1; D. armeniaca 3_04 Da2 + second run; D. unisexualis 224 Du1); in each read both
windows were located (megablast, best piece >= 150 bp at >= 85% reaching the window end next to the copy); the allele is
the read sequence between the projected copy ends, CLsat in it by ssearch36 (>= 60 bits); each read was assigned to the
assembly whose place region it matches best (hybrids: subgenome). Toy control: the place regions themselves as reads
returned the assembly answers for all 8 places of Sval00049 and Smix00074.
`reads.tsv`: animal, id, read, strand, window identities, allele_bp, CLsat_bp, best assembly, best and second bits.
