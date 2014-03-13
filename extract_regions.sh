#!/bin/bash
# Generate appropriate bed files to use with the metaplotting...
#
master=/Users/stephanie/ll/results/FP/dREG_regions_marked.bed.gz
folder=/Users/stephanie/ll/results/FP

gunzip -c $master > m.temp

# region-divisions
awk 'BEGIN{OFS="\t"}{ if ($13==1) print $1,$2,$3,$4,$16,$18 }' m.temp | sed '1d' > $folder/dREG_gs.bed
# restrict gene body hits to be >2k from the gene start, for the plotting
awk 'BEGIN{OFS="\t"}{ if (($14!=0)&&($17>2000)) print $1,$2,$3,$4,$16,$18 }' m.temp | sed '1d' > $folder/dREG_gb.bed
# non-genes don't have strands!
awk 'BEGIN{OFS="\t"}{ if ($15==1) print $1,$2,$3,$4,$16 }' m.temp | sed '1d' > $folder/dREG_ng.bed
