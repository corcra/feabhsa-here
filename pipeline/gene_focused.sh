#!/bin/bash

suffix=40m_SVM      # eg 40m_SVM, full_SVM
#suffix=full_SVM      # eg 40m_SVM, full_SVM

THRE_high=0.8
THRE_low=0.5

genefolder=/Users/stephanie/ll/data/genes

results=/Users/stephanie/ll/results/$suffix
confident=$results/dREG_regions_confident_$THRE_high.bed.gz
maybe=$results/dREG_regions_maybe_$THRE_low.bed.gz

genelist=$genefolder/gene_list.bed

# --- get the gene-dREG overlap for both confidence datasets... --- #
# # --- this may be unnecessary based on shifting this part of the analysis to get_covariates... let's see... -#
gunzip -c $maybe | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - > $results/genes_maybe_overlap.bed
python gene_focused_tidy.py $results/genes_maybe_overlap.bed $results/genes_maybe_overlap_tidy.bed
mv $results/genes_maybe_overlap_tidy.bed $results/genes_maybe_overlap.bed

gunzip -c $confident | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - > $results/genes_confident_overlap.bed
python gene_focused_tidy.py $results/genes_confident_overlap.bed $results/genes_confident_overlap_tidy.bed
mv $results/genes_confident_overlap_tidy.bed $results/genes_confident_overlap.bed

# --- produce a list of genes with hits in their TSS - based on the maybe data --- #
awk '{ if ($8!=0) print $0 }' $results/genes_maybe_overlap.bed > $results/active_genes.bed

# --- produce a list of genes with hits in their gene body - based on the confident data --- #
awk '{ if ($9!=0) print $0 }' $results/genes_confident_overlap.bed | bedmap --skip-unmapped --echo - $results/active_genes.bed > $results/active_genes_gb_hit.bed

# --- produce a list of genes with NO HITS in their gene body - based on the maybe data --- #
awk '{ if ($9==0) print $0 }' $results/active_genes.bed > $results/active_genes_nogb_hit.bed

echo "# genes with dREG hit at TSS:" `wc -l $results/active_genes.bed | awk '{ print $1 }'`
# now, which of these have hits in the body?
echo "# genes with dREG hit at TSS and in gene body:" `wc -l $results/active_genes_gb_hit.bed | awk '{ print $1 }'`

# --- create a list for running the HMM... --- #
cat $results/active_genes_gb_hit.bed $results/active_genes_nogb_hit.bed > $results/genes_for_analysis.bed

# --- clean up... don't really need these after the master genelist has been created... --- #
#rm -v active_genes.bed
#rm -v genes_gs_gb_hit.bed
#rm -v genes_gs_nogb_hit.bed

## --- produce the goodlist --- #
#gunzip -c $confident | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - > $results/gene_list_with_overlaps.bed
#python gene_focused_tidy.py $results/gene_list_with_overlaps.bed $results/gene_list_with_overlaps_tidy.bed
#mv $results/gene_list_with_overlaps_tidy.bed $results/gene_list_with_overlaps.bed
# --- produce the bad list --- #
#gunzip -c $maybe | sed '1d' | bedops --not-element-of -1
#

# get the genes with hits at gene_start
#awk '{ if ($7==0) print $0 }' $results/gene_list_with_overlaps.bed > $results/genes_no_hits.bed


# not sure how useful thi sstuff is #
# get the genes with no hits at gene_start
#awk '{ if (($7>0)&&($8==0)) print $0 }' $results/gene_list_with_overlaps.bed > $results/genes_only_body_hits.bed
#awk '{ print $10 }' $results/genes_only_body_hits.bed | tr ';' '\n' > $results/dREG_IDs_onlybody.txt
#echo -e 'ID\tstart\tbody' > $results/gene_IDs_bodyhit.txt
#awk 'BEGIN{OFS="\t"}{ if ($9>0) print $4,$8,$9 }' $results/gene_list_with_overlaps.bed >> $results/gene_IDs_bodyhit.txt
# (use this in the visualisation script! when are these guys discovered?)


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
