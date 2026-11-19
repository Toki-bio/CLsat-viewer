#!/bin/bash
# Reads behind the parental block test, for viewing (owner 2026-11-19). Same reads and mapping as checks2.sh part B
# (first 30 M R1 reads of D. raddei nairensis 17_78 and D. valentini 5_38, minimap2 -x sr, unn = unp.fna), but kept as
# alignments in three windows of unn contig JAWDKN010000070.1: nairensis-type sequence (18.50 Mb), valentini-type
# sequence after the assembly switch (19.50 Mb), and the start of the 334-kb array L17 (22,777,927).
# usage: blockview.sh <toy|full>
set -euo pipefail
W=/data/V/toki/clsat_map/naicheck
cd $W || exit 1
export TMPDIR=$W/tmp; mkdir -p tmp
O=blockview_$1; mkdir -p $O
R=/data/V/toki/Darevskia_reads; MM=/data/V/toki/bin/minimap2
if [ "$1" = toy ]; then NB=20000; else NB=30000000; fi
C=JAWDKN010000070.1
printf '%s\t18495000\t18505000\tnai_part\n%s\t19495000\t19505000\tval_part\n%s\t22770000\t22790000\tL17_start\n' $C $C $C > $O/windows.bed
for s in nai_17_78:$R/raddei_nairensis/2014/D-raddei-nairensis_17_78_1.fastq.gz val_5_38:$R/valentini/2014/D-valentini_5_38_1.fastq.gz; do
  n=${s%%:*}; f=${s#*:}
  set +o pipefail
  zcat $f | head -n $((NB*4)) | $MM -ax sr --secondary=no -t 32 -R "@RG\tID:$n\tSM:$n" idx/unp.sr.mmi - 2>/dev/null \
    | samtools view -u -L $O/windows.bed - | samtools sort -o $O/$n.unn.bam -
  set -o pipefail
  samtools index $O/$n.unn.bam
  echo "$n $(samtools view -c $O/$n.unn.bam) reads in windows"
done
samtools faidx /data/V/toki/Darevskia_v2/genomes/unp.fna $C > $O/$C.fa && samtools faidx $O/$C.fa   # whole contig, 23.1 Mb, so BAM coordinates match
echo ALLDONE
