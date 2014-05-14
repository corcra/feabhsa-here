#!/bin/bash

# The purpose of this is to - given a list of genes with associated transition point info (thus timepoint), etc.

suffix=40m_SVM
THRE_high=0.8
THRE_low=0.5
promoter_size=1000
time_region=$1
version=C

# --- datasets yo --- #
covar_folder=/Users/stephanie/ll/results/rate_analysis/$time_region
genes=/Users/stephanie/ll/results/$suffix/active_genes.bed
nicegenes=${genes/.bed/_Rable.bed}
awk '{ print $1, $2, $3, $4, $5, $6 }' $genes > $nicegenes
confident=/Users/stephanie/ll/results/$suffix/dREG_regions_confident_$THRE_high.bed.gz
maybe=/Users/stephanie/ll/results/$suffix/dREG_regions_maybe_$THRE_low.bed.gz
introns=/Users/stephanie/ll/data/genes/intronic.bed
exons=/Users/stephanie/ll/data/genes/exonic.bed
# --- note, using the plus strand for the CpG data, but it doesn't matter... just gives positive values
CpG=/Users/stephanie/ll/data/ChIPseq/CpG_plus.bedgraph
H3K79me2=/Users/stephanie/ll/data/ChIPseq/H3K79me2_sorted.bed.gz

if [ $time_region == "early" ]
then
    from=/Users/stephanie/ll/from_hojoong/result/$version\_5min_bin500.40m_SVM.txt
    to=/Users/stephanie/ll/from_hojoong/result/$version\_1230min_bin1000.40m_SVM.txt
fi

if [ $time_region == "mid" ]
then
    from=/Users/stephanie/ll/from_hojoong/result/$version\_1230min_bin1000.40m_SVM.txt
    to=/Users/stephanie/ll/from_hojoong/result/$version\_25min_bin2000.40m_SVM.txt
fi

if [ $time_region == "late" ]
then
    from=/Users/stephanie/ll/from_hojoong/result/$version\_25min_bin2000.40m_SVM.txt
    to=/Users/stephanie/ll/from_hojoong/result/$version\_50min_bin5000.40m_SVM.txt
fi

# --- create the 'region of relevance' --- #
R --slave --file=get_FP_affected_region.r --args $nicegenes $from $to
sort-bed regions.bed.temp > region_of_relevance.bed
#awk 'BEGIN{OFS="\t"}{ if ($6=="+") print $1, $2+'"$from"',$2+'"$to"',$4,$5,$6; else  print $1, $3-'"$to"', $3-'"$from"',$4,$5,$6  }' $genes  > region_of_relevance.bed

# --- number of dREG hits in region of relevance --- #
echo "name" "n_gb" > $covar_folder/n_gb.txt
#gunzip -c $confident | sed '1d' | bedmap --range 250 --echo --count --skip-unmapped --delim "\t" region_of_relevance.bed - | awk '{ print $4,$7}' >> $covar_folder/n_gb.txt
gunzip -c $confident | sed '1d' | bedmap --range 250 --echo --count --delim "\t" region_of_relevance.bed - | awk '{ print $4,$6*1000/($3-$2)}' >> $covar_folder/n_gb.txt
## if there are no hits with the maybe dataset, count it as no hits...
#gunzip -c $maybe | sed '1d' | bedmap --range 250 --echo --indicator region_of_relevance.bed - | grep '|0' | awk '{ print $4,0}' >> $covar_folder/n_gb.txt

# --- length of intron 1 --- #
# note, this is produced by process_refseq.py
# note, this is for all genes - need to specify to specific gene, etc.
echo "name" "int1len" > $covar_folder/int1len.txt
awk '{ if ($5==1) print $4, $3-$2 }' $introns >> $covar_folder/int1len.txt

# --- number of exons in region of relevance --- #
echo "name" "n_exon" > $covar_folder/n_exon.txt
bedmap --echo --echo-map-id --delim ";" region_of_relevance.bed $exons | python get_exon_counts.py >> $covar_folder/n_exon.txt

# --- nucleosome density at promoter region --- # ... or is it just position of first nucleosome?
echo "name" "nuc_pos" > $covar_folder/nuc_pos.txt
# for now, ...?

# --- CpG content in the region --- #
echo "name" "CpG" > $covar_folder/CpG_gb.txt
bedmap --echo --sum --delim "\t" region_of_relevance.bed $CpG | awk '{ print $4, $6*1000/($3-$2) }' >> $covar_folder/CpG_gb.txt

# --- H3K79me2 --- #
echo "name" "H3K79me2" > $covar_folder/H3K79me2_gb.txt
gunzip -c $H3K79me2 | bedmap --echo --count --delim "\t" region_of_relevance.bed - | awk '{ print $4, $6*1000/($3-$2) }' >> $covar_folder/H3K79me2_gb.txt

# -- tidy up --- #
rm *.temp
