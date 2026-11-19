# Parental block test of the Darevskia hybrid assemblies (CLsat revision, 2026-11-19)

`index.html` - the page (opens at https://toki-bio.github.io/CLsat-viewer/block_test/): every 100-kb block of unv, unn,
armf and armm with the differences of two parental-species animals' short reads; blocks called for the other parent
and the one CLsat array affected (unn L17); the unn contig JAWDKN010000070.1 block by block (assembly switch at 18.8 Mb);
reads in three windows of that contig.

Rebuild: `python build.py` (writes short_read_blocks.html; index.html = that file wrapped in html/body tags). Inputs:
- `B_blocks_called.tsv`, `B_<animal>_on_<unm|unp>.tsv` - D. unisexualis block tables (KIT /data/V/toki/clsat_map/naicheck/out2_full/; unm = unv, unp = unn)
- `arm/called_*.tsv` - D. armeniaca block calls, two pairs of animals (KIT /data/V/toki/clsat_map/armsplit/full/; column 4 = D. mixta reads, 5 = D. valentini reads)
- `reads.tsv` - per-read difference positions in the three windows (`blockview.sh`, `blockview_reads.sh`, run on KIT; BAMs in KIT /data/V/toki/clsat_map/naicheck/blockview_full/)
- locus table: CLsat_workflow run_gb_v2 regions/locus_master.tsv (path set in build.py)

Rule (Methods §1): reads with mapping quality >= 30; a block is called for the animal whose reads differ less, margin >= 25%, >= 300 reads of each animal.
