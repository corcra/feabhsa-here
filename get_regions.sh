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

# Find out where the hits are!
echo "Dividing hits into regions."
bedmap --fraction-ref 0.5 --indicator pred.temp $TSS > TSS.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb.temp
bedmap --fraction-ref 0.5 --indicator pred.temp $non_gene > non_gene.temp

# Get the distance to the nearest TSS... (for all hits!) (modify the master)
wc -l pred.temp
wc -l gene_starts.temp
bedtools closest -D "a" -a pred.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-1), $NF }' > dist.pre.temp
wc -l dist.temp

echo "Combining!"
paste pred.temp TSS.temp gb.temp non_gene.temp dist.temp | gzip -c > $predfile.ehhh

echo "Tidying up!"
#rm -v *.temp
