#!/bin/bash
# author:       Stephanie Hyland (sh985@cornell.edu)
# date:         February 2014
# description:  Take a list of dREG (chr start end) and assign to genomic regions. Also, get distance to closest gene, and histone marks.

data=/Users/stephanie/ll/data/

# Genes
gene_list=$data/genes/gene_list.bed
gene_start=$data/genes/gene_start.bed
gene_body=$data/genes/gene_body.bed
non_gene=$data/genes/non_genes.bed
awk 'BEGIN{OFS="\t"}{ if ($6=="+") { print $1, $2, $2+1, $4, $6 } else { print $1, $3-1, $3, $4, $6 } }' $gene_list > gene_starts.temp

# Markers (note using bedGraphs)
H3K4me1=$data/ChIPseq/H3K4me1_Meiss.sorted.bed.gz
H3K4me3=$data/ChIPseq/H3K4me3_Marson.sorted.bed.gz
H3K27ac=$data/ChIPseq/H3K27Ac.sorted.bed.gz

# You know what, just give me the whole path.
predfile=$1
gunzip -c $predfile | sed '/inflation/d' > pred.temp
# (stripping the header line)

# Find out where the hits are!
echo "Dividing hits into regions."
bedmap --fraction-ref 0.5 --indicator pred.temp $gene_start > gene_start.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb.temp
bedmap --fraction-ref 0.5 --indicator pred.temp $non_gene > non_gene.temp

# Get the distance to the nearest gene_start... (for all hits!) (modify the master)
echo "Getting distance to closest gene start."
bedtools closest -d -a pred.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-2), $NF, $(NF-1) }' > dist.pre.temp
python post_bedtools_join_closest.py
# (creates dist.temp)

# Get histone counts! (totally unnormalised here...)
echo "Getting histone counts."
gunzip -c $H3K4me1 | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K4me1.temp
gunzip -c $H3K4me3 | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K4me3.temp
gunzip -c $H3K27ac | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K27ac.temp

echo "Combining!"
paste pred.temp gene_start.temp gb.temp non_gene.temp dist.temp H3K4me1.temp H3K4me3.temp H3K27ac.temp > comb.temp
gunzip -c $predfile | head -1 | awk '{ print $0, "gene_start","gene_body","non_gene","closest_gene","distance","strand","H3K4me1","H3K4me3","H3K27ac" }' > headerfile.temp
cat headerfile.temp comb.temp | gzip -c > ${predfile/_uniq/_marked}

echo "Tidying up!"
#rm -v *.temp
