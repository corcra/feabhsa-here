#!/bin/bash
# Preparing the gene annotations.
# Note, for the purpose of finding TSS, we will extend the gene window back by 500, first.

cd ~/ll/data/genes/

# Create genes from alt4.txt - use NM name as ID!
awk 'BEGIN{OFS="\t"}{ print $3, $5, $6, $2, "NA", $4 }' alt4.txt | sort-bed - > gene_list_with_dupes.bed
echo "Removing duplicates (same start and end points, same chromosome, same strand)."
python ~/Projects/feabhsa-here/remove_duplicate_isoforms.py gene_list_with_dupes.bed

# Extend
echo "Extending to include TSS."
awk 'BEGIN{OFS="\t"}{ if ($6=="+") print $1,$2-500,$3,$4,$5,$6; else print $1,$2,$3+501,$4,$5,$6 }' gene_list.bed | sort-bed - > genes_inc_TSS.bed

# Get the negative set - non-gene regions (also excludes TSS)
# This just means looking for anywhere in the mm9.bed that isn't in genes_inc_TSS.bed
bedops --difference mm9.bed genes_inc_TSS.bed | awk '{ print $1, $2, $3, "non_gene" }' > non_genes.bed

# Look at strand-specific non-gene regions (not using this right now).
# Split into positive and negative strands.
echo "Splitting strands."
sed -n '/+/p' genes_inc_TSS.bed > genes_inc_TSS_plus.bed
sed -n '/-/p' genes_inc_TSS.bed > genes_inc_TSS_minus.bed
echo "Getting negative set."
bedops --difference mm9.bed genes_inc_TSS_plus.bed | awk '{ print $1, $2, $3, "non_gene", "NA", "+" }' > non_genes_plus.bed 
bedops --difference mm9.bed genes_inc_TSS_minus.bed | awk '{ print $1, $2, $3, "non_gene", "NA", "-" }' > non_genes_minus.bed 
cat non_genes_plus.bed non_genes_minus.bed | sort-bed - > non_genes_with_strand.bed

# To make genes_TSS
echo "Getting TSS."
awk 'BEGIN{OFS="\t"}{ if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$3-500,$3+501,$4,$5,$6 }' gene_list.bed > TSS.bed

# Genes_no_TSS
echo "Getting gene body."
awk 'BEGIN{OFS="\t"}{ if ($6=="+") print $1,$2+501,$3,$4,$5,$6; else print $1,$2,$3-500,$4,$5,$6 }' gene_list.bed | awk '{ if (($3-$2)>500) print $0 }' | sort-bed - > gene_body_with_TSS.bed
echo "Subtracting TSS inside genes."
bedmap --delim "|" --multidelim "|" --echo --echo-map gene_body_with_TSS.bed TSS.bed > gbt.temp
python ~/Projects/feabhsa-here/remove_TSS.py gbt.temp 
# (produces gene_body.bed)

# TSS inside genes
echo "Looking for TSS hits inside gene bodies."
bedmap --delim "|" --multidelim "|" --echo --echo-map --skip-unmapped --bases-uniq-f TSS.bed gene_body_with_TSS.bed > TSS_in_gene_body.bed

# Record results.
echo "Recording results."
echo `date` > logfile.txt
echo "Total gene list:" `wc -l gene_list_with_dupes.bed` >> logfile.txt
echo "Unique gene list:" `wc -l gene_list.bed` >> logfile.txt
echo "TSS list:" `wc -l TSS.bed` >> logfile.txt
echo "Gene body:" `wc -l gene_body.bed` >> logfile.txt
echo "Non gene regions (+):" `wc -l non_genes_plus.bed` >> logfile.txt
echo "Non gene regions (-):" `wc -l non_genes_minus.bed` >> logfile.txt
echo "TSS with gene body overlap (either strand):" `wc -l TSS_in_gene_body.bed` >> logfile.txt
echo "TSS entirely contained:" `awk 'BEGIN{FS="|"}{ if($NF==1) print $0 }' TSS_in_gene_body.bed | wc -l` >> logfile.txt

echo "Tidying!"
rm -v gene_list_with_dupes.bed
rm -v genes_inc_TSS*.bed
rm -v non_genes_minus.bed
rm -v non_genes_plus.bed
rm -v gbt.temp
rm -v gene_body_with_TSS.bed
