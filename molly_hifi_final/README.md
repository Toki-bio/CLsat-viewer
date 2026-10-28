# Amazon molly satellite families, rerun with the final workflow (stages 1-8)

Same four HiFi assemblies as before (lat = *P. latipinna*, mex = *P. mexicana*, Flat / Fmex = Amazon molly copies inherited from latipinna / mexicana), run through satpipe (CLsat workflow) at the final workflow version, `discovery: {strand: both, monomerize: true}`. 102 cross-species clusters, 69 with at least 3 members. Cluster names are new (c000...), they do not match the old famNNN numbers. Each cluster name carries its member counts: `c010_n62_len187_lat..._Fmex...` = cluster 10, 62 members, 187 bp consensus, members per genome.

Groups are the author's manual peeling of the old run, mapped onto the new clusters by sequence (best ssearch36 -z 3 hit). Consensus = majority of a MAFFT alignment of the cluster members (members are locus-level monomer consensuses). Start points are NOT unified between groups: different groups are different families, so there is no common frame; check start points by eye.

- [All 69 family consensuses aligned together](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fall_families_aligned.fa&title=Molly%20final%20all%20family%20consensuses)

| group | what | clusters | ViewAlign |
|---|---|---|---|
| G1 | 186 bp group (old 013,014,015,017,019,028) | c010, c011, c017 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FG1.fa&title=Molly%20final%20G1) |
| G2 | 83-84 bp group (old 067-120, 28 families) | c052, c053, c057, c058, c059, c060, c061, c062, c063, c064, c065, c066, c067, c068, c069, c070, c071, c072, c073, c074, c076, c077, c080, c082, c084, c086, c088 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FG2.fa&title=Molly%20final%20G2) |
| G3 | 151-152 bp (old 047,048,049) | c028 | (one cluster) |
| G4 | 105-111 bp, less certain (old 060,063,064) | c038, c039 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FG4.fa&title=Molly%20final%20G4) |
| G5 | 80 bp singleton (old 135) | c098 | (one cluster) |
| G6 | 250-256 bp pair (old 001,002) | c000 | (one cluster) |
| G7 | 186 bp pair (old 016,018) | c012 | (one cluster) |
| G8 | 166 bp pair (old 040,041) | c026 | (one cluster) |
| G9 | 176 bp group (old 030,032,033,034) | c018, c019 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FG9.fa&title=Molly%20final%20G9) |
| G10 | 183 bp pair (old 020,024) | c013 | (one cluster) |
| G11 | 175 bp singleton (old 035) | c016 | (one cluster) |
| G12 | 100 bp singleton (old 062) | c040 | (one cluster) |
| G13 | 94 bp singleton (old 045) | c095 | (one cluster) |
| G14 | 80-94 bp pair (old 051,133) | c101 | (one cluster) |
| G15 | 83 bp singleton (old 107) | c055 | (one cluster) |
| G16a | 175 bp singleton (old 031) | c020 | (one cluster) |
| G16b | 171-175 bp pair (old 036,039) | c021, c025 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FG16b.fa&title=Molly%20final%20G16b) |
| G17 | 151 bp singleton (old 050) | c029 | (one cluster) |
| G18 | 122 bp singleton (old 055) | c034 | (one cluster) |
| G19 | 110 bp singleton (old 059) | c037 | (one cluster) |
| G20 | 112 bp singleton (old 058) | c036 | (one cluster) |
| GX | new clusters with >=3 members not in any of the author's old groups (to check) | c002, c009, c014, c015, c022, c027, c033, c041, c042, c043, c044, c046, c051, c056, c081, c089, c094 | [open](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fgroups%2FGX.fa&title=Molly%20final%20GX) |

Relations between length classes (all-vs-all, both strands, ssearch36 -z 3): no significant hit between the 84, 105, 151, 176, 186 and 250 bp classes; 310 strong within-class pairs, only 1 reverse complement.

Files: `all_families_aligned.fa`, `family_consensus.fa`, `groups/*.fa`.

## MG02 vs MG15: one family or two? (2026-10-28)

[Open the alignment in ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Fmolly_hifi_final%2Fmg02_mg15.aln.fa&title=Molly%20MG02%20vs%20MG15) (129 rows, `mafft --adjustdirection`; `_R_` = reverse complemented).
Rows, top to bottom by name prefix: `MG02_consensus_84bp`, `MG15_consensus_82bp`; `MG02_clus_*` / `MG15_clus_*` = the run's cluster consensuses behind each family (MG15 = c055 only, old fam107); then 3 monomers per locus:
`BOTH_*` = 12 MG02 loci where MG15 also scores >= 50 bits; `MG15_*` = all 9 MG15 loci; `MG02_*` = 12 random MG02 loci where MG15 scores below half of MG02.
Row name = class_locus_genome_bestFamBits_secondFamBits_monomer. Locus list: [selected.tsv](selected.tsv).
MG15 is the runner-up at all 305 MG02 loci and MG02 at all 9 MG15 loci. The two consensuses align only over 64 bp at 78 % identity.
