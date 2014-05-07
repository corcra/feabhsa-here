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
gunzip -c $maybe | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - | python gene_focused_tidy.py > $results/genes_maybe_overlap.bed

gunzip -c $confident | sed '1d' | bedmap --range 500 --multidelim "|" --delim "|" --echo --count --echo-map $genelist - | python gene_focused_tidy.py > $results/genes_confident_overlap.bed

# --- produce a list of genes with hits in their TSS - based on the maybe data --- #
awk '{ if ($8!=0) print $0 }' $results/genes_maybe_overlap.bed > $results/active_genes.bed

# PROBS DON'T NEED ANY OF THIS #
# --- produce a list of genes with hits in their gene body - based on the confident data --- #
awk '{ if ($9!=0) print $0 }' $results/genes_confident_overlap.bed | bedmap --skip-unmapped --echo - $results/active_genes.bed > $results/active_genes_gb_hit.bed

# --- produce a list of genes with NO HITS in their gene body - based on the maybe data --- #
awk '{ if ($9==0) print $0 }' $results/active_genes.bed > $results/active_genes_nogb_hit.bed

echo "# genes with dREG hit at TSS:" `wc -l $results/active_genes.bed | awk '{ print $1 }'`
# now, which of these have hits in the body?
echo "# genes with dREG hit at TSS and in gene body:" `wc -l $results/active_genes_gb_hit.bed | awk '{ print $1 }'`
# --- create a list for running the HMM... --- #
#cat $results/active_genes_gb_hit.bed $results/active_genes_nogb_hit.bed > $results/genes_for_analysis.bed
