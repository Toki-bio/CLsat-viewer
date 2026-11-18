#!/bin/bash
# Trace dispersed CLsat solo copies to the same place in the other five assemblies, by single-copy flanks.
# Input: dispersed.tsv (id asm contig start end family subfamily; 1-based). Genomes: genomes/<asm>.fa (+ .fai).
# Rules: flank = the nearest 400-bp window (0-4.6 kb from the copy, 200-bp steps) that is single-copy in its own assembly
# (one hit >= 200 bp at >= 90%); its distance is added back when the inner end is projected in the target;
# in a target, the flank's inner end must be in the best hit (second best < 0.9 x best bits); a place = both flanks on
# one contig, same orientation, in order, inner ends <= 50 kb apart. Between the inner ends: CLsat by ssearch36 against
# the doubled family consensus (>= 60 bits).
set -euo pipefail
W=/data/V/toki/clsat_map/solo_orth_2026-11-18
cd "$W" || exit 1
IN=${1:-dispersed.tsv}; OUT=${2:-out}; T=${THREADS:-24}
WIN=400; SHOW=300
ASMS="val mix unv unn armf armm"
mkdir -p db "$OUT"/{fl,hit,reg,aln}
for a in $ASMS; do [ -s db/$a.nsq ] || [ -s db/$a.00.nsq ] || makeblastdb -in genomes/$a.fa -dbtype nucl -out db/$a > /dev/null; done

rc() { rev | tr ACGTacgtNn TGCAtgcaNn; }
seqof() { samtools faidx genomes/$1.fa "$2:$3-$4" | grep -v '>' | tr -d '\n'; }

: > "$OUT/status.tsv"; : > "$OUT/regions.tsv"; : > "$OUT/flanks_used.tsv"
while IFS=$'\t' read -r id a c s e fam sub; do
  L=$(awk -v c="$c" '$1==c{print $2}' genomes/$a.fa.fai)
  # candidate flank windows: WIN bp at distance d = 0, 200, ... , 4600 from the copy on each side
  : > "$OUT/fl/$id.cand.fa"
  for d in $(seq 0 200 4600); do
    l2=$((s-1-d)); l1=$((l2-WIN+1)); r1=$((e+1+d)); r2=$((r1+WIN-1))
    [ $l1 -ge 1 ] && { echo ">L_$d"; seqof $a "$c" $l1 $l2; echo; } >> "$OUT/fl/$id.cand.fa"
    [ $r2 -le $L ] && { echo ">R_$d"; seqof $a "$c" $r1 $r2; echo; } >> "$OUT/fl/$id.cand.fa"
  done
  # single-copy test in its own assembly: exactly one hit >= half the window at >= 90%; keep the nearest such window
  blastn -task megablast -query "$OUT/fl/$id.cand.fa" -db db/$a -perc_identity 90 -evalue 1e-20 -num_threads $T \
         -outfmt '6 qseqid length' -max_target_seqs 50 | awk -v h=$((WIN/2)) '$2>=h{n[$1]++} END{for(k in n) print k"\t"n[k]}' \
         > "$OUT/fl/$id.cn.tsv"
  dL=$(awk -F'\t' '$1~/^L_/ && $2==1{split($1,x,"_"); print x[2]}' "$OUT/fl/$id.cn.tsv" | sort -n | head -1)
  dR=$(awk -F'\t' '$1~/^R_/ && $2==1{split($1,x,"_"); print x[2]}' "$OUT/fl/$id.cn.tsv" | sort -n | head -1)
  okL=0; okR=0; : > "$OUT/fl/$id.fa"
  [ -n "$dL" ] && { okL=1; echo ">${id}_L"; seqof $a "$c" $((s-dL-WIN)) $((s-1-dL)); echo; } >> "$OUT/fl/$id.fa"
  [ -n "$dR" ] && { okR=1; echo ">${id}_R"; seqof $a "$c" $((e+1+dR)) $((e+dR+WIN)); echo; } >> "$OUT/fl/$id.fa"
  dL=${dL:-0}; dR=${dR:-0}; qL=$WIN
  printf '%s\tL_window_distance=%s\tR_window_distance=%s\tL_unique=%s\tR_unique=%s\n' "$id" "$dL" "$dR" "$okL" "$okR" >> "$OUT/flanks_used.tsv"
  for t in $ASMS; do
    [ "$t" = "$a" ] && continue
    # per flank: best hit by bits; accept if it contains the inner end (within 50 bp) and second best < 0.9 x best
    [ -s "$OUT/fl/$id.fa" ] || { : > "$OUT/hit/$id.$t.tsv"; }
    [ -s "$OUT/fl/$id.fa" ] && blastn -task megablast -query "$OUT/fl/$id.fa" -db db/$t -evalue 1e-20 -num_threads $T -max_target_seqs 20 \
      -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore' > "$OUT/hit/$id.$t.tsv"
    res=$(awk -v oL=$okL -v oR=$okR -v qL=$qL -v dL=$dL -v dR=$dR -v id=$id '
      { k=$1; if(!(k in b1) || $9>b1[k]){ b2[k]=b1[k]; b1[k]=$9; h[k]=$0 } else if($9>b2[k]) b2[k]=$9 }
      END{
        for(side=0; side<2; side++){
          k = side==0 ? id"_L" : id"_R"; ok = side==0 ? oL : oR
          if(!ok){ out[side]="nonunique"; continue }
          if(!(k in h)){ out[side]="nohit"; continue }
          if(b2[k]!="" && b2[k]>=0.9*b1[k]){ out[side]="ambiguous"; continue }
          split(h[k],x,"\t"); qs=x[5]; qe=x[6]; ss=x[7]; se=x[8]; str=(se>=ss)?"+":"-"
          # inner end projected to the copy boundary over the window distance (dL, dR)
          if(side==0){ if(qe<qL-50){ out[side]="innerend_missing"; continue }; tin = (str=="+") ? se+(qL-qe)+dL : se-(qL-qe)-dL }
          else       { if(qs>51){ out[side]="innerend_missing"; continue };    tin = (str=="+") ? ss-(qs-1)-dR : ss+(qs-1)+dR }
          out[side]=x[2]":"tin":"str":"x[3]
        }
        print out[0]"\t"out[1]
      }' "$OUT/hit/$id.$t.tsv")
    hL=$(echo "$res" | cut -f1); hR=$(echo "$res" | cut -f2)
    status=""; gap=NA; clsat=0; reg=""
    if [[ "$hL" == *:*:* && "$hR" == *:*:* ]]; then
      IFS=: read -r cL pL sL iL <<< "$hL"; IFS=: read -r cR pR sR iR <<< "$hR"
      if [ "$cL" = "$cR" ] && [ "$sL" = "$sR" ]; then
        if [ "$sL" = "+" ]; then gap=$((pR-pL-1)); a1=$((pL+1)); a2=$((pR-1)); else gap=$((pL-pR-1)); a1=$((pR+1)); a2=$((pL-1)); fi
        if [ $gap -ge -50 ] && [ $gap -le 50000 ]; then
          if [ $gap -gt 0 ]; then
            { echo ">gap"; seqof $t "$cL" $a1 $a2; echo; } > "$OUT/reg/$id.$t.gap.fa"
            clsat=$(ssearch36 -n -z 3 -E 1000 -m 8C -T 4 q/${fam}x2.fa "$OUT/reg/$id.$t.gap.fa" 2>/dev/null | \
                    awk '!/^#/ && $12>=60 {s=($9<$10)?$9:$10; e=($9<$10)?$10:$9; print s"\t"e}' | sort -k1,1n | \
                    awk '{ if($1>ce){ tot+=ce-cs; cs=$1; ce=$2 } else if($2>ce) ce=$2 } END{ tot+=ce-cs; print tot+0 }')
          fi
          if [ "$clsat" -gt 0 ]; then status=present; elif [ $gap -le 50 ]; then status=empty; else status=other_sequence; fi
          # region for the alignment: SHOW bp of each flank + the gap, in the source orientation
          if [ "$sL" = "+" ]; then r1=$((pL-SHOW+1)); r2=$((pR+SHOW-1)); else r1=$((pR-SHOW+1)); r2=$((pL+SHOW-1)); fi
          reg="$cL:$r1-$r2:$sL"
        else status=flanks_too_far; fi
      else status=flanks_disagree; fi
    elif [[ "$hL" == *:*:* || "$hR" == *:*:* ]]; then status=one_flank_only
      # region for the alignment: SHOW bp of the found flank + copy length + SHOW bp beyond, in the source orientation
      n=$((e-s+1))
      if [[ "$hL" == *:*:* ]]; then IFS=: read -r cX pX sX iX <<< "$hL"
        if [ "$sX" = "+" ]; then r1=$((pX-SHOW+1)); r2=$((pX+n+SHOW)); else r1=$((pX-n-SHOW)); r2=$((pX+SHOW-1)); fi
      else IFS=: read -r cX pX sX iX <<< "$hR"
        if [ "$sX" = "+" ]; then r1=$((pX-n-SHOW)); r2=$((pX+SHOW-1)); else r1=$((pX-SHOW+1)); r2=$((pX+n+SHOW)); fi
      fi
      [ $r1 -lt 1 ] && r1=1
      reg="$cX:$r1-$r2:$sX:one_flank"
    else status=not_found; fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$a" "$sub" "$((e-s+1))" "$t" "$hL" "$hR" "$status" "$gap" "$clsat" >> "$OUT/status.tsv"
    [ -n "$reg" ] && printf '%s\t%s\t%s\n' "$id" "$t" "$reg" >> "$OUT/regions.tsv"
  done
  # alignment: source copy with SHOW bp flanks, then every target place found by both flanks
  { echo ">${id}|${a}|${c}:$((s-SHOW))-$((e+SHOW))|${sub}|source"; seqof $a "$c" $((s-SHOW)) $((e+SHOW)); echo
    awk -v id=$id '$1==id' "$OUT/regions.tsv" 2>/dev/null | while IFS=$'\t' read -r _ t r; do
      IFS=: read -r rc_ rr st kind <<< "$r"; r1=${rr%-*}; r2=${rr#*-}
      stt=$(awk -F'\t' -v id=$id -v t=$t '$1==id && $5==t{print $8 (($10>0)?"_CLsat"$10"bp":"")}' "$OUT/status.tsv")
      echo ">${id}|${t}|${rc_}:${r1}-${r2}|${st}|${stt}"
      if [ "$st" = "+" ]; then seqof $t "$rc_" $r1 $r2; else seqof $t "$rc_" $r1 $r2 | rc; fi; echo
    done; } > "$OUT/aln/$id.in.fa"
  n=$(grep -c '>' "$OUT/aln/$id.in.fa")
  if [ "$n" -ge 2 ]; then mafft --quiet --thread 8 --localpair --maxiterate 1000 "$OUT/aln/$id.in.fa" > "$OUT/aln/$id.aln.fa"; fi
  echo "$id done ($n sequences)" >&2
done < "$IN"
