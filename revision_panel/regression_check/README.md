# Published-panel regression: disputed loci

Regression of the final workflow (fixed ssearch36 statistics, all audit fixes) on the six GenBank assemblies of the
submitted manuscript (KIT run 2026-10-25). 527 of 627 segmented loci are identical; these are the loci whose family
changed or whose split differs. Each alignment: the four family consensuses, then consecutive one-unit (146 bp)
windows of the region (each a rotation of one unit), turned to the CLsat1 frame, named
`contig:position|best family|identity|bits` (ssearch36 -z 3 against the doubled consensuses).

| region | published | final workflow | windows | windows support | alignment |
|---|---|---|---|---|---|
| D. armeniaca JAWWNF010000515.1:197000-237029 | CLsat3 | CLsat1 | 200 of 200 windows CLsat1 | new | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam1_arm_JAWWNF010000515_197000-237029.aln.fa&title=D.%20armeniaca%20JAWWNF010000515.1%3A197000-237029) |
| D. mixta JAWWNH010000118.1:40000-50925 | CLsat1 | CLsat2 | 73 of 74 windows CLsat1 | published | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam2_mix_JAWWNH010000118_40000-50925.aln.fa&title=D.%20mixta%20JAWWNH010000118.1%3A40000-50925) |
| D. mixta JAWWNH010000265.1:0-10732 | CLsat1 | CLsat2 | 52 CLsat2, 15 CLsat1, 6 CLsat4 | new | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam3_mix_JAWWNH010000265_0-10732.aln.fa&title=D.%20mixta%20JAWWNH010000265.1%3A0-10732) |
| D. valentini JAWWNG010000030.1:355000-360655 (published small CLsat3 locus L494, now inside a CLsat1 locus) | CLsat3 | CLsat1 | 19 CLsat3, 10 CLsat2, 6 CLsat1, 3 CLsat4 | unclear | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam4_val_JAWWNG010000030_355000-360655.aln.fa&title=D.%20valentini%20JAWWNG010000030.1%3A355000-360655%20%28published%20small%20CLsat3%20locus%20L494%2C%20now%20inside%20a%20CLsat1%20locus%29) |
| D. valentini JAWWNG010000555.1:91000-653030 (562 kb, 200 of ~3850 windows) | CLsat3 | CLsat1 | 197 of 200 windows CLsat1 | new | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam5_val_JAWWNG010000555_91000-653030.aln.fa&title=D.%20valentini%20JAWWNG010000555.1%3A91000-653030%20%28562%20kb%2C%20200%20of%20~3850%20windows%29) |
| D. valentini JAWWNG010003007.1:38000-132028 | CLsat1 | CLsat3 | 199 of 200 windows CLsat1 | published | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Ffam6_val_JAWWNG010003007_38000-132028.aln.fa&title=D.%20valentini%20JAWWNG010003007.1%3A38000-132028) |
| D. armeniaca JAWWNF010000280.1:1005000-1018119 (published end 1018119, new end 1010773) |  |  | 51 CLsat3, 16 CLsat2, 13 CLsat1, 9 CLsat4 | see alignment | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Fend_arm_JAWWNF010000280_1005000-1018119.aln.fa&title=D.%20armeniaca%20JAWWNF010000280.1%3A1005000-1018119%20%28published%20end%201018119%2C%20new%20end%201010773%29) |
| unp JAWDKN010000113.1:75000-130000 (published split at 85000 and 120000, new one locus) |  |  | 81 CLsat1, 45 CLsat4, 42 CLsat2, 32 CLsat3 | see alignment | [ViewAlign](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FCLsat-viewer%2Fmain%2Frevision_panel%2Fregression_check%2Fsplit_unp_JAWDKN010000113_75000-130000.aln.fa&title=unp%20JAWDKN010000113.1%3A75000-130000%20%28published%20split%20at%2085000%20and%20120000%2C%20new%20one%20locus%29) |

Arbitration assigns a segment the family with the highest single window score; one strong window of another family
can decide a whole segment (fam2, fam6 here; fam5 in the published run, by 0.4 bits). Segmentation itself uses the
family that wins each 1-kb bin. In the final-panel catalogue 11 of 554 loci (0.66 Mb) get a family that differs
from the one winning most of their bins.
