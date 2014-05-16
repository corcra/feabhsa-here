#!/bin/bash
echo "get_covar_minimal.sh: reporting for duty!"

# --- options --- #
suffix=40m_SVM
THRE_high=0.8
THRE_low=0.5
promoter_size=1000
time_region=$1

# --- paths --- #
covar_folder=/Users/stephanie/ll/results/rate_analysis/$time_region

# --- datasets yo --- #
confident=/Users/stephanie/ll/results/$suffix/dREG_regions_confident_$THRE_high.bed.gz
maybe=/Users/stephanie/ll/results/$suffix/dREG_regions_maybe_$THRE_low.bed.gz
region=/Users/stephanie/Projects/feabhsa-here/pipeline/$time_region\_region_of_relevance.bed

# --- raw covars --- # (note, using plus for CpG, but it's all the same)
introns=/Users/stephanie/ll/data/genes/intronic.bed
exons=/Users/stephanie/ll/data/genes/exonic.bed
CpG=/Users/stephanie/ll/data/ChIPseq/CpG_plus.bedgraph
H3K79me2=/Users/stephanie/ll/data/ChIPseq/H3K79me2_sorted.bed.gz

# --- number of dREG hits in region of relevance --- #
echo "get_covar_minimal.sh: n_gb!"
echo "name" "n_gb" "d_gb"> $covar_folder/n_gb.txt
#gunzip -c $confident | sed '1d' | bedmap --range 250 --echo --count --skip-unmapped --delim "\t" $region - | awk '{ print $4,$7}' >> $covar_folder/n_gb.txt
gunzip -c $confident | sed '1d' | bedmap --range 250 --echo --count --delim "\t" $region - | awk '{ print $4,$6,$6*1000/($3-$2)}' >> $covar_folder/n_gb.txt
## if there are no hits with the maybe dataset, count it as no hits...
#gunzip -c $maybe | sed '1d' | bedmap --range 250 --echo --indicator $region - | grep '|0' | awk '{ print $4,0}' >> $covar_folder/n_gb.txt

# --- length of intron 1 --- #
echo "get_covar_minimal.sh: int1len!"
# note, this is produced by process_refseq.py
# note, this is for all genes - need to specify to specific gene, etc.
echo "name" "int1len" > $covar_folder/int1len.txt
awk '{ if ($5==1) print $4, $3-$2 }' $introns >> $covar_folder/int1len.txt

# --- number of exons in region of relevance --- #
echo "get_covar_minimal.sh: n_exon!"
echo "name" "n_exon" "d_exon"> $covar_folder/n_exon.txt
bedmap --echo --echo-map-id --delim ";" $region $exons | python get_exon_counts.py >> $covar_folder/n_exon.txt

# --- nucleosome density at promoter region --- # ... or is it just position of first nucleosome?
echo "get_covar_minimal.sh: nuc_pos!"
echo "name" "nuc_pos" > $covar_folder/nuc_pos.txt
# for now, ...?

# --- CpG content in the region --- #
echo "get_covar_minimal.sh: CpG!"
echo "name" "CpG" > $covar_folder/CpG_gb.txt
bedmap --echo --sum --delim "\t" $region $CpG | awk '{ print $4, $6*1000/($3-$2) }' >> $covar_folder/CpG_gb.txt

# --- H3K79me2 --- #
echo "get_covar_minimal.sh: H3K79me2!"
echo "name" "H3K79me2" > $covar_folder/H3K79me2_gb.txt
gunzip -c $H3K79me2 | bedmap --echo --count --delim "\t" $region - | awk '{ print $4, $6*1000/($3-$2) }' >> $covar_folder/H3K79me2_gb.txt

# -- tidy up --- #
rm *.temp
