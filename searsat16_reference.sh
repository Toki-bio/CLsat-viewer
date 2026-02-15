#!/bin/bash
# searsat16 - Split a satellite locus sequence into blocks of tandem repeats
# Usage: searsat16 <query.fa> <locus.fa> [Slf] [%Homology]
# $1 - query consensus
# $2 - bank (locus sequence)
# $3 - Slf (sequence length fraction, default 0.9)
# $4 - %Homology (default 75)
#
# Key features:
# - Uses ssearch36 Smith-Waterman to find tandem repeat matches
# - Splits locus into separate tandem arrays based on gaps (>100bp)
# - Only processes arrays with ≥4 tandems
# - Aligns each array independently with MAFFT
#
# ver.14: added ssearch36 -r +2/-2 parameter to improve hanging ends of monomer sequences

if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied!'
    exit 1
fi

Files=($(which "ssearch36") $(which "bedtools") $(which "seqkit") $(which "seqret") "$1" )
for f in "${Files[@]}" ; do
    if [ ! -f "$f" ]; then
        echo "File $f: not found"
        exit 1
    fi
done

if [ -z "$3" ]
    then Slf=0.9
    else Slf=$3
fi

len=$(awk '/^>/{if (l!="") print l; l=0; next}{l+=length($0)}END{print l}' $1) # Query length calculation
flen=$(awk -v a=$Slf -v b=$len 'BEGIN{if (a>1) flen=a; else flen=int(a*b); print flen;}') # Overlap length calculation
bname=$2

if [ -z "$4" ]
    then Homology=75
    else Homology=$4
fi

# FIRST and the only ROUND
echo Searching $1 in $2
ssearch36 -r +2/-2 -g -3 -T 24 -Q -n -N 5000 -z 11 -E 2 -w 95 -W 70 -m 8C $1 $2 3 \
    | awk  '/^[^#]/ {printf($2"\t"$9-1"\t"$10"\t"$10-$9+1"\t"); if ($8-$7>0) print$8-$7+1"\t""+""\t"$3; else print-($8-$7+1)"\t""-""\t"$3 }' \
    | sort -V \
    | bedtools merge -s -d -1 -c 4,5,6,7 -o max,max,distinct,max \
    | awk '{if ($7>='$Homology' && $4>='$flen') print}' \
    | sort -V > plus_$2-$1.bed \
    || exit 1

awk '{a[$6]++;} END{for(i in a) print a[i]"  "i}' plus_$2-$1.bed\
    | awk -v bname=$bname '{if ($2=="-") {print bname,$1" reverse complement hits found"}}' >> $2.ssat.log
awk '$6!="-"' plus_$2-$1.bed > $2.hits.bed # filtering out rc hits

bl=$(awk 'END {print NR}' $2.hits.bed) # counting the number of all hits in a locus

if [[ "$bl" -lt "3" ]] # checking if more than 3 "+" hits were found
    then
        echo Skipping
        echo "$2 $bl plus hits in a locus with less than 3 tandems found, skipped" >> $2.ssat.log
        rm plus_$2-$1.bed
        rm $2.hits.bed
        exit 111
    else
        echo "Splitting locus"
fi

# CRITICAL STEP: Split hits into separate arrays based on 100bp gap threshold
cat $2.hits.bed | awk 'OFS="\t" {xa[NR]=$2}{xb[NR]=$3} END {for (i in xa) print $1,xa[i],xb[i],xa[i]-xb[i-1],$4,$5,$6}' \
    | awk 'OFS="\t" {if ($4>100 && NR>1) print"\n"$1,$2,$3,$5,$6,$7; else print $1,$2,$3,$5,$6,$7}' \
    | awk -v RS= -v bname=$bname '{print > (bname".hits-" NR ".bed")}' # Splitting hits into separate bed files

# Analyze each group of possible tandems in a locus
for b in $bname.hits-*.bed; do
    bn=$b
    wc -l <$b | awk -v bn=$bn -v bname=$bname '{print bname,bn, $0, "hits found in a piece of locus"}' >> $2.ssat.log
    CheckFile=$(awk 'END {print NR, FILENAME}' $b | awk '{if ($1<4) print $2; else print "1"}')
    
    if [[ "$CheckFile" == "$b" ]] # Skip if <4 hits in array
        then
            echo "Less than 4 hits found in" $b
        else
            echo "Analysing hits"
            bnl=$(awk 'END {print NR}' $b)
            blen=$(awk '/^>/{if (l!="") print l; l=0; next}{l+=length($0)}END{print l}' $2)
            
            # Recalculate coordinates to make continuous sequence, expand boundaries ±100bp
            tac $b \
                | awk 'OFS="\t" {xa[NR]=$2}{xb[NR]=$3} END {for (i in xa) print $1,xa[i],xb[i],xa[i-1]-1,$4,$5,$6}' \
                | awk 'OFS="\t" {if (NR==1) print $1,$2,$3,$5,$6,$7; else print $1,$2,$4,$5,$6,$7}' \
                | awk -v blen=$blen -v bnl=$bnl 'OFS="\t" {if (NR==1 && $3+99<blen) $3=$3+100}{if (NR==1 && $3+99>=blen) $3=blen}{if (NR==bnl && $2-99>0) $2=$2-100}{if (NR==bnl && $2-99<=0) $2=0} 1' \
                | awk 'OFS="\t" {if (NR>1) $3=$3+1; print $0}' | tac > h_$b.bed
            
            # Extract sequences
            bedtools getfasta -s -fi $2 -bed h_$b.bed > $2_$b.bnk
            
            awk '{printf $0}' $b | awk '{print $1,$2,$(NF-3),$(NF-2),$(NF-1),$(NF)}' | awk '{gsub(/:|-|\(/," ")} {print}' | awk -v bnl=$bnl 'OFS="\t" {if ($4==")") print $1,$3-$6,$3-$5,$7,$8,"-",bnl; else print $1,$2+$5,$2+$6,$7,$8,"+",bnl}' >> $2.matches.bed
            
            rm h_$b.bed
            sed -i -e '$a\' $1
            
            # Combine consensus + extracted sequences, align with MAFFT
            cat $1 $2_$b.bnk > $2_$b_c.bnk
            mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --quiet $2_$b_c.bnk | seqkit seq -w 0 > $bnl$b$2.al
            rm $2_$b_c.bnk
        fi
done

awk '/bed|plus/ {if (NF==10) print $3,$1; else {print $2,$1}}' $2.ssat.log > $2.splus.log
sort -k 1,1n $2.splus.log | awk '{a[$1]++;} END{for(i in a) print a[i]"  "i}' > $2.count.log
rm $2*hits*.bed plus_$2* $2_* $2*matches*.bed $2*.fai $2*.log

echo "Done"
