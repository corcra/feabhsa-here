#!/bin/bash

gunzip -c dREG_regions_marked.bed.gz | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map gene_list.bed - > gene_list_with_overlaps.bed
python gene_focused_tidy.py

# R code
# library(ggplot2)
# dat<-read.table('gene_list_with_overlaps_tidy.bed')
# names(dat)<-c("chr","start","end","gene","score","strand","total","promoter","body")
# ggplot(dat,aes(x=total))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in genes.")
# ggsave("hits_ingene.pdf")
# ggplot(dat,aes(x=promoter))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in promoters.")
# ggsave("hits_inpromoter.pdf")
# ggplot(dat,aes(x=body))+geom_histogram(binwidth=0.5)+theme_bw()+ggtitle("Number of dREG hits in gene bodies.")
# ggsave("hits_inbody.pdf")
#

# get the genes with no hits at TSS
awk '{ if ($7==0) print $0 }' gene_list_with_overlaps_tidy.bed > genes_no_hits.bed
awk '{ if (($7>0)&&($8==0)) print $0 }' gene_list_with_overlaps_tidy.bed > genes_only_body_hits.bed
awk '{ print $10 }' genes_only_body_hits.bed | tr ';' '\n' > dREG_IDs_onlybody.txt
# (use this in the visualisation script! when are these guys discovered?)
