#!/bin/bash

master=/Users/stephanie/ll/results/FP/dREG_regions_marked.bed.gz
genefolder=/Users/stephanie/ll/data/genes
geneanalysis=$genefolder/analyse_FP
genelist=$genefolder/gene_list.bed
gunzip -c $master | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - > $geneanalysis/gene_list_with_overlaps.bed
python gene_focused_tidy.py
mv $geneanalysis/gene_list_with_overlaps_tidy.bed $geneanalysis/gene_list_with_overlaps.bed

# R code
# library(ggplot2)
# dat<-read.table('gene_list_with_overlaps_tidy.bed')
# names(dat)<-c("chr","start","end","gene","score","strand","total","promoter","body")
# ggplot(dat,aes(x=total))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in genes.")
# ggsave("pdfs/hits_ingene.pdf")
# ggplot(dat,aes(x=promoter))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in promoters.")
# ggsave("pdfs/hits_inpromoter.pdf")
# ggplot(dat,aes(x=body))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in gene bodies.")
# ggsave("pdfs/hits_inbody.pdf")
#

# get the genes with no hits at gene_start
awk '{ if ($7==0) print $0 }' $geneanalysis/gene_list_with_overlaps.bed > $geneanalysis/genes_no_hits.bed
awk '{ if (($7>0)&&($8==0)) print $0 }' $geneanalysis/gene_list_with_overlaps.bed > $geneanalysis/genes_only_body_hits.bed
awk '{ print $10 }' $geneanalysis/genes_only_body_hits.bed | tr ';' '\n' > $geneanalysis/dREG_IDs_onlybody.txt
echo -e 'ID\tstart\tbody' > $geneanalysis/gene_IDs_bodyhit.txt
awk 'BEGIN{OFS="\t"}{ if ($9>0) print $4,$8,$9 }' $geneanalysis/gene_list_with_overlaps.bed >> $geneanalysis/gene_IDs_bodyhit.txt
# (use this in the visualisation script! when are these guys discovered?)
