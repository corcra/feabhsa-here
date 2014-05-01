#!/bin/bash

# The purpose of this is to - given a list of genes with associated transition point info (thus timepoint), etc.

suffix=40m_SVM
THRE_high=0.8
THRE_low=0.5

covar_folder=/Users/stephanie/ll/results/rate_analysis
genes=/Users/stephanie/ll/results/$suffix/genes_for_analysis.bed
confident=/Users/stephanie/ll/results/$suffix/dREG_regions_confident_$THRE_high.bed.gz
maybe=/Users/stephanie/ll/results/$suffix/dREG_regions_maybe_$THRE_low.bed.gz

time_region=$1

if [ $time_region == "early" ]
then
    from=0
    to=7500
fi

if [ $time_region == "mid" ]
then
    from=7500
    to=30000
fi

if [ $time_region == "late" ]
then
    from=30000
    to=100000
fi

# --- create the 'region of relevance' --- #
awk 'BEGIN{OFS="\t"}{ print $1, $2+'"$from"',$2+'"$to"',$4,$5,$6 }' $genes  > region_of_relevance.bed

# --- number of dREG hits in region of relevance --- #
echo "name" "n_gb" > $covar_folder/n_gb.txt
gunzip -c $confident | sed '1d' | bedmap --range 500 --echo --count --skip-unmapped --delim "\t" region_of_relevance.bed - | awk '{ print $4,$7}' >> $covar_folder/n_gb.txt
# if there are no hits with the maybe dataset, count it as no hits...
gunzip -c $maybe | sed '1d' | bedmap --range 500 --echo --indicator region_of_relevance.bed - | grep '|0' | awk '{ print $4,0}' >> $covar_folder/n_gb.txt

# --- length of intron 1 --- #
# note, this is produced by process_refseq.py
# note, this is for all genes - need to specify to specific gene, etc.
introns=/Users/stephanie/ll/data/genes/intronic.bed
echo "name" "int1len" > $covar_folder/int1len.txt
awk '{ if ($5==1) print $4, $3-$2 }' $introns >> $covar_folder/int1len.txt

# --- number of exons in region of relevance --- #
exons=/Users/stephanie/ll/data/genes/exonic.bed
echo "name" "n_exon" > $covar_folder/n_exon.txt
bedmap --echo --count --delim "\t" region_of_relevance.bed $exons | awk '{ print $4, $7 }' >> $covar_folder/n_exon.txt
