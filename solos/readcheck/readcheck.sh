#!/bin/bash
# Read check of the traced solo places with each animal's own HiFi reads (no assembly used in the test itself).
# For each place: the two single-copy 400-bp flank windows chosen by trace.sh (out_all/fl/<id>.fa, distances in
# out_all/flanks_used.tsv). Reads holding window words are fished; in each read both windows are located (megablast,
# best piece >= 150 bp at >= 85%, reaching the window end next to the copy); the allele = sequence between the projected copy ends; CLsat in it by ssearch36 (>= 60 bits).
# In a hybrid, each spanning read is assigned to the assembly whose place region it matches best (aln/<id>.in.fa).
# Usage: readcheck.sh IDS_FILE OUTDIR "animal:bam[,bam]" ...     (toy: animal:fasta)
set -euo pipefail
W=/data/V/toki/clsat_map/solo_orth_2026-11-18
cd "$W" || exit 1
IDS=$1; OUT=$2; shift 2
mkdir -p "$OUT"
TR=out_all
# 1. words: 25-mers every 8 bp from both windows of every place, both strands, no low-complexity
: > "$OUT/windows.fa"
while read -r id; do cat "$TR/fl/$id.fa" >> "$OUT/windows.fa"; done < "$IDS"
awk '/^>/{next}{s=toupper($0); for(i=1;i+24<=length(s);i+=8){w=substr(s,i,25); if(w~/N/) continue;
      n["A"]=gsub(/A/,"A",w); n["C"]=gsub(/C/,"C",w); n["G"]=gsub(/G/,"G",w); n["T"]=gsub(/T/,"T",w);
      if(n["A"]>14||n["C"]>14||n["G"]>14||n["T"]>14) continue; if(w~/AAAAAA|CCCCCC|GGGGGG|TTTTTT/) continue; print w}}' \
  "$OUT/windows.fa" | sort -u > "$OUT/words1.txt"
rev "$OUT/words1.txt" | tr ACGT TGCA | cat "$OUT/words1.txt" - | sort -u > "$OUT/words.txt"
echo "words: $(wc -l < "$OUT/words.txt")" >&2
# 2. fish per animal (in parallel): reads holding >= 1 word; spanning is decided by blast below
for spec in "$@"; do
  animal=${spec%%:*}; files=${spec#*:}
  (
    for f in ${files//,/ }; do
      case "$f" in
        *.bam) samtools view -@ 4 "$f" | awk -F'\t' '{print ">"$1"\n"$10}' ;;
        *)     awk '/^>/{if(s)print h"\n"s; h=$1; s=""; next}{s=s $0} END{if(s)print h"\n"s}' "$f" ;;
      esac
    done | grep -B1 --no-group-separator -F -f "$OUT/words.txt" > "$OUT/reads_$animal.fa" || true
    echo "$animal fished $(grep -c '>' "$OUT/reads_$animal.fa")" >&2
  ) &
done
wait
# 3. per animal and place: spanning reads, allele length, CLsat, assignment
printf 'animal\tid\tread\tstrand\twinL_id\twinR_id\tallele_bp\tCLsat_bp\tbest_assembly\tbest_bits\tsecond_bits\n' > "$OUT/reads.tsv"
for spec in "$@"; do
  animal=${spec%%:*}; R="$OUT/reads_$animal.fa"; [ -s "$R" ] || continue
  makeblastdb -in "$R" -dbtype nucl -out "$OUT/db_$animal" > /dev/null; samtools faidx "$R"
  while read -r id; do
    fam=$(awk -F'\t' -v id=$id '$1==id{print $6}' dispersed.tsv)
    dL=$(awk -F'\t' -v id=$id '$1==id{split($2,x,"="); print x[2]}' $TR/flanks_used.tsv)
    dR=$(awk -F'\t' -v id=$id '$1==id{split($3,x,"="); print x[2]}' $TR/flanks_used.tsv)
    blastn -task megablast -query "$TR/fl/$id.fa" -db "$OUT/db_$animal" -evalue 1e-20 -perc_identity 85 -max_target_seqs 5000 \
      -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore slen' 2>/dev/null | \
    awk -v id=$id -v dL=$dL -v dR=$dR '$4>=150 {
        k=$2 SUBSEP $1; if(!(k in b) || $9>b[k]){ b[k]=$9; h[k]=$0 } }
      END{ for(k in h){ split(k,kk,SUBSEP); r=kk[1]; side=(kk[2]~/_L$/)?"L":"R"; H[r,side]=h[k]; R_[r]=1 }
        for(r in R_){ if(!((r,"L") in H) || !((r,"R") in H)) continue
          split(H[r,"L"],a,"\t"); split(H[r,"R"],c,"\t")
          sa=(a[8]>=a[7])?"+":"-"; sc=(c[8]>=c[7])?"+":"-"; if(sa!=sc) continue
          if(a[6]<350 || c[5]>51) continue      # best piece must reach the window end next to the copy
          # inner ends projected over the window distance; L window end = q 400, R window start = q 1
          if(sa=="+"){ pL=a[8]+(400-a[6])+dL; pR=c[7]-(c[5]-1)-dR; al=pR-pL-1; a1=pL+1; a2=pR-1 }
          else       { pL=a[8]-(400-a[6])-dL; pR=c[7]+(c[5]-1)+dR; al=pL-pR-1; a1=pR+1; a2=pL-1 }
          if(al<-100 || al>60000) continue
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n", id, r, sa, a[3], c[3], al, a1, a2 } }' > "$OUT/span.$animal.$id.tsv"
    while IFS=$'\t' read -r _ r st iL iR al a1 a2; do
      cl=0
      if [ "$al" -gt 20 ]; then
        samtools faidx "$R" "$r:$a1-$a2" > "$OUT/tmp.$animal.fa"
        cl=$(ssearch36 -n -z 3 -E 1000 -m 8C -T 4 q/${fam}x2.fa "$OUT/tmp.$animal.fa" 2>/dev/null | \
             awk '!/^#/ && $12>=60 {s=($9<$10)?$9:$10; e=($9<$10)?$10:$9; print s"\t"e}' | sort -k1,1n | \
             awk '{ if($1>ce){ tot+=ce-cs; cs=$1; ce=$2 } else if($2>ce) ce=$2 } END{ tot+=ce-cs; print tot+0 }')
      fi
      # assignment: read segment (windows + allele) against every place region of this copy (aln/<id>.in.fa)
      lo=$(( a1-500 < 1 ? 1 : a1-500 )); hi=$(( a2+500 ))
      samtools faidx "$R" "$r:$lo-$hi" > "$OUT/seg.$animal.fa" 2>/dev/null
      asg=$(blastn -task megablast -query "$OUT/seg.$animal.fa" -subject "$TR/aln/$id.in.fa" -outfmt '6 sseqid bitscore' 2>/dev/null | \
            awk '{split($1,x,"|"); a=x[2]; s[a]+=$2} END{for(k in s) print s[k]"\t"k}' | sort -k1,1gr | head -2 | \
            awk 'NR==1{b=$2; bb=$1} NR==2{sb=$1} END{printf "%s\t%s\t%s", b, bb, sb+0}')
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$animal" "$id" "$r" "$st" "$iL" "$iR" "$al" "$cl" "$asg" >> "$OUT/reads.tsv"
    done < "$OUT/span.$animal.$id.tsv"
  done < "$IDS"
done
echo done >&2
