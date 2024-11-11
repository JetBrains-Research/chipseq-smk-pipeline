#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }

>&2 echo "bam2tagsbw.sh $@"
if [[ $# -lt 4 ]]; then
    echo "Need 4 parameters! <UNIQUE_BAM> <FRAGMENT> <CHROM_SIZES> <OUTPUT_BW>"
    exit 1
fi

UNIQUE_BAM=$1
FRAGMENT=$2
SHIFT=$(($FRAGMENT / 2))
CHROM_SIZES=$3
OUTPUT_BW=$4

T=$'\t'

bedtools bamtobed -i ${UNIQUE_BAM} |\
  grep -E "chr[0-9XYM]+$T" |\
  awk -v OFS='\t' -v S=${SHIFT} \
    '{if ($6 != "-") {print($1, $2+S, $2+S+1)} else {if ($3-S>=1) {print($1, $3-S, $3-S+1)}}}' |\
    sort -u -k1,1 -k3,3n -k2,2n > ${OUTPUT_BW}.bed_fragment;

cat ${OUTPUT_BW}.bed_fragment |\
  awk -v OFS='\t' 'BEGIN{C="";S=0;E=0;X=0}
    {if(C!=$1||S!=$2||E!=$3){if(X!=0){print(C,S,E,X)};C=$1;S=$2;E=$3;X=1}else{X=X+1}}
    END{if(X!=0){print(C,S,E,X)}}' > ${OUTPUT_BW}.bdg;
rm ${OUTPUT_BW}.bed_fragment;

# Remove coordinates outside chromosome sizes
bedtools slop -i ${OUTPUT_BW}.bdg -g ${CHROM_SIZES} -b 0 | bedClip stdin ${CHROM_SIZES} ${OUTPUT_BW}.bdg.clip
rm ${OUTPUT_BW}.bdg;

# Fix problem with not sorted clip file
LC_COLLATE=C sort -k1,1 -k2,2n ${OUTPUT_BW}.bdg.clip > ${OUTPUT_BW}.bdg.clip.sorted;
rm ${OUTPUT_BW}.bdg.clip;

bedGraphToBigWig ${OUTPUT_BW}.bdg.clip.sorted ${CHROM_SIZES} ${OUTPUT_BW};
rm ${OUTPUT_BW}.bdg.clip.sorted;