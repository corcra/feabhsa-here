#!/bin/bash
# author:       Stephanie Hyland (sh985@cornell.edu)
# date:         February 2014
# description:  Take a list of dREG (chr start end) and assign to genomic regions. Also, get distance to closest gene.

# Gene data
data=/Users/stephanie/ll/data/
gene_list=$data/genes/gene_list.bed
TSS=$data/genes/TSS.bed
gene_body=$data/genes/gene_body.bed
non_gene=$data/genes/non_genes.bed
awk 'BEGIN{OFS="\t"}{ print $1, $2, $2+1, $4 }' $gene_list > gene_starts.temp

# You know what, just give me the whole path.
predfile=$1
gunzip -c $predfile | sed '/inflation/d' > pred.temp
# (stripping the header line)

# Find out where the hits are!
echo "Dividing hits into regions."
bedmap --fraction-ref 0.5 --indicator pred.temp $TSS > TSS.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb.temp
bedmap --fraction-ref 0.5 --indicator pred.temp $non_gene > non_gene.temp

# Get the distance to the nearest TSS... (for all hits!) (modify the master)
bedtools closest -D "a" -a pred.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-1), $NF }' > dist.pre.temp
python post_bedtools_join_closest.py
# (creates dist.temp)

echo "Combining!"
paste pred.temp TSS.temp gb.temp non_gene.temp dist.temp > comb.temp
gunzip -c $predfile | head -1 | awk '{ print $0, "TSS","gene_body","non_gene","closest_gene","distance" }' > headerfile.temp
cat headerfile.temp comb.temp | gzip -c > $predfile.ehhh

# Header line

echo "Tidying up!"
rm -v *.temp
