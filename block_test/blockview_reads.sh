#!/bin/bash
# Per-read difference positions in the three windows of blockview.sh (for the page): window, animal, start, end (1-based
# reference), mapping quality, reference positions of mismatches, deletions and insertions (comma list).
set -euo pipefail
W=/data/V/toki/clsat_map/naicheck
cd $W/blockview_full || exit 1
C=JAWDKN010000070.1
: > reads.tsv
while read -r c s e name; do
  for an in nai_17_78 val_5_38; do
    samtools view -h $an.unn.bam "$c:$((s+1))-$e" | samtools calmd - $C.fa 2>/dev/null | \
    awk -v w=$name -v a=${an%%_*} -v s=$s -v e=$e 'BEGIN{OFS="\t"} /^@/{next} {
      if(and($2,4)) next
      pos=$4; cig=$6; md=""; for(i=12;i<=NF;i++) if($i~/^MD:Z:/) md=substr($i,6)
      # reference length and insertion positions from CIGAR
      rl=0; ins=""; p=pos
      while(match(cig,/^[0-9]+[MIDNSHP=X]/)){n=substr(cig,1,RLENGTH-1)+0; op=substr(cig,RLENGTH,1); cig=substr(cig,RLENGTH+1)
        if(op=="M"||op=="D"||op=="N"||op=="="||op=="X"){rl+=n; p+=n} else if(op=="I") ins=ins (ins==""?"":",") p }
      end=pos+rl-1; if(end<=s || pos>e) next
      # mismatches and deletions from MD
      out=ins; p=pos
      while(md!=""){
        if(match(md,/^[0-9]+/)){p+=substr(md,1,RLENGTH)+0; md=substr(md,RLENGTH+1)}
        else if(substr(md,1,1)=="^"){md=substr(md,2); match(md,/^[A-Z]+/); out=out (out==""?"":",") p; p+=RLENGTH; md=substr(md,RLENGTH+1)}
        else {out=out (out==""?"":",") p; p++; md=substr(md,2)}
      }
      print w, a, pos, end, $5, out }' >> reads.tsv
  done
done < windows.bed
awk -F'\t' '{n[$1" "$2]++; q[$1" "$2]+=($5>=30)} END{for(k in n) print k, n[k], "reads,", q[k], "with MAPQ>=30"}' reads.tsv | sort
