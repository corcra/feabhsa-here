#!/bin/bash

suffix=40m_SVM      # eg 40m_SVM, full_SVM
#suffix=full_SVM      # eg 40m_SVM, full_SVM

genefolder=/Users/stephanie/ll/data/genes

results=/Users/stephanie/ll/results/$suffix
confident=$results/FP_dREG_confident_0.8.bed.gz
uncertain=$results/FP_dREG_uncertain_0.5.bed.gz
genelist=$genefolder/gene_list.bed

# --- produce the goodlist --- #
gunzip -c $confident | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - > $results/gene_list_with_overlaps.bed
python gene_focused_tidy.py $results/gene_list_with_overlaps.bed $results/gene_list_with_overlaps_tidy.bed
mv $results/gene_list_with_overlaps_tidy.bed $results/gene_list_with_overlaps.bed
# --- produce the bad list --- #
#gunzip -c $uncertain | sed '1d' | bedops --not-element-of -1

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

# get the genes with hits at gene_start
awk '{ if ($7!=0) print $0 }' $results/gene_list_with_overlaps.bed > $results/genes_gs_hit.bed
echo "# genes with dREG hit at TSS:" `wc -l $results/genes_gs_hit.bed | awk '{ print $1 }'`
# now, which of these have hits in the body?
echo "# genes with dREG hit at TSS and in gene body:" `awk '{ if ($8!=0) print $0 }' $results/genes_gs_hit.bed | wc -l`


# not sure how useful thi sstuff is #
# get the genes with no hits at gene_start
awk '{ if ($7==0) print $0 }' $results/gene_list_with_overlaps.bed > $results/genes_no_hits.bed
awk '{ if (($7>0)&&($8==0)) print $0 }' $results/gene_list_with_overlaps.bed > $results/genes_only_body_hits.bed
awk '{ print $10 }' $results/genes_only_body_hits.bed | tr ';' '\n' > $results/dREG_IDs_onlybody.txt
echo -e 'ID\tstart\tbody' > $results/gene_IDs_bodyhit.txt
awk 'BEGIN{OFS="\t"}{ if ($9>0) print $4,$8,$9 }' $results/gene_list_with_overlaps.bed >> $results/gene_IDs_bodyhit.txt
# (use this in the visualisation script! when are these guys discovered?)
