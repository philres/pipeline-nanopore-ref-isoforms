#!/bin/bash - 
set -o nounset                              # Treat unset variables as an error

CTX=20
REF=$1
CORES=$2
IN_BAM=$3
OUT_BAM=$4
TMP="filter_internal_$$"
TMP="filter_internal_139323"

mkdir -p $TMP
seqkit bam -j $CORES -f Read,Ref,Pos,EndPos,Strand $IN_BAM 2> $TMP/aln_info.tsv
csvtk -t -j $CORES mutate2 -L 0 -n Start -e "\$Strand == 1 ? \$EndPos - $CTX : \$Pos - $CTX" $TMP/aln_info.tsv |\
csvtk -t -j $CORES mutate2 -L 0 -n End -e "\$Strand == 1 ? \$EndPos + $CTX : \$Pos + $CTX" - |\
csvtk -t -j $CORES mutate2 -L 0 -n BedStrand -e "\$Strand == 1 ? '+' : '-' " - |\
csvtk -t -j $CORES mutate2 -L 0 -n BedScore -e "1" - |\
csvtk -t -j $CORES cut -f Ref,Start,End,Read,BedScore,BedStrand - | tail -n +2 > $TMP/three_prime.bed
bedtools getfasta -nameOnly -s -bed $TMP/three_prime.bed -fi $REF -tab |\
sed 's/([-\+])\t/\t/'> $TMP/three_prime.tsv
awk '{ tmp=$2;gsub(/[TGC]+/,"", tmp); if (length(tmp) > 20 ) {print $1}}' $TMP/three_prime.tsv > $TMP/bad_reads.tsv
seqkit -j $CORES bam -x -Q -f Strand -G $TMP/bad_reads.tsv $IN_BAM 2>/dev/null > $OUT_BAM
echo "Filtered out " `wc -l $TMP/bad_reads.tsv | cut -d $' ' -f 1` "internal priming reads."
rm -r $TMP
