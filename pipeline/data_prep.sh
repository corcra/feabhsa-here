#!/bin/bash
# author:       Stephanie Hyland (sh985@cornell.edu)
# date:         February-April 2014
# description:  This is something of a 'master' pipeline. It does the following:
#               1. Combine dREG predictions (above thresh) from all timepoints (into $raw_high).
#               2. Format this nicely (via process_master.py)
#               3. Classify all points by genomic region, get overlap with histone marks.

# Threshold ... obviously, select these carefully
THRE_high=0.8
THRE_low=0.5

# Inputs!
datatype=FP     # eg FP, TRP... not sure why I'm still asking this
suffix=40m_SVM  # eg full_SVM, 40m_SVM... other?

# Paths
data=/Users/stephanie/ll/data
folder=$data/dREG/$suffix
results=/Users/stephanie/ll/results/$suffix

# Main files (note - no QC on the maybe hits! only use this to exclude things, really...)
# --- confident ---
raw_high=$results/dREG_regions_$THRE_high.bed.gz
proc_high=$results/dREG_regions_proc_$THRE_high.bed.gz
pre_QC_confident=$results/dREG_regions_preQC_$THRE_high.bed.gz
confident=$results/dREG_regions_confident_$THRE_high.bed.gz
# --- uncertain ---
raw_low=$results/dREG_regions_$THRE_low.bed.gz
proc_low=$results/dREG_regions_maybe_proc_$THRE_low.bed.gz
maybe=$results/dREG_regions_maybe_$THRE_low.bed.gz

# Genes (for later)
gene_list=$data/genes/gene_list.bed
gene_start=$data/genes/gene_start.bed
gene_body=$data/genes/gene_body.bed
non_gene=$data/genes/non_genes.bed

# Markers (also for later)
H3K4me1=$data/ChIPseq/H3K4me1_Meiss.sorted.bed.gz
H3K4me3=$data/ChIPseq/H3K4me3_Marson.sorted.bed.gz
H3K27ac=$data/ChIPseq/H3K27Ac.sorted.bed.gz

# Mappability
mappable=$data/mappability/mm9map.bedgraph

# --- Let's do this! --- #
for datatime in 0min 2min 5min 12.5min 25min 50min
do
    echo "Encorporating" $datatime
    # Set up vars
    preds=$folder/$datatype\_$datatime.predictions.bedGraph.gz

    # remove multiply-mapped regions ... actually, should be doing this pre-dREG!
    # using map-id seems weird, but in the bedgraph the value goes in the fourth column
    gunzip -c $preds | bedmap --skip-unmapped --echo --echo-map-id - $mappable | awk 'BEGIN{FS="|"}{ if ($2==1) print $1 }' > $preds.temp
    gunzip -c $preds > $preds.temp
    high_pos=$results/$datatype\_$datatime.pos.$THRE_high.bedGraph.gz
    low_pos=$results/$datatype\_$datatime.pos.$THRE_low.bedGraph.gz

    # What is the score distribution? (can plot this later as necessary)
    if ! [ -a ${preds/.predictions/.gs.predictions} ]
    then
        echo "Getting scores."
        bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $gene_start | gzip -c > ${preds/.predictions/.gs.predictions}
        bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $gene_body | gzip -c > ${preds/.predictions/.gb.predictions}
        bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $non_gene | gzip -c > ${preds/.predictions/.ng.predictions}
    fi

    # Expands positive regions into window +-50, merges these.
    echo "Getting high-confidence regions."
    awk 'BEGIN{OFS="\t"} ($4 > '"$THRE_high"') { print $1,$2-50,$3+51,$4 }' $preds.temp | sort-bed - | bedops --merge - | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, "'$datatype'_" "'${datatime/min/}'_" NR }' > $high_pos.temp
    echo "Getting low-confidence regions."
    awk 'BEGIN{OFS="\t"} ($4 > '"$THRE_low"') { print $1,$2-50,$3+51,$4 }' $preds.temp | sort-bed - | bedops --merge - | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, "'$datatype'_" "'${datatime/min/}'_" NR }' > $low_pos.temp

    rm -v $preds.temp

    # --- high-confidence --- #
    # Create the master list, if necessary
    if ! [ -a $raw_high ]
    then
        echo "Creating master list (high conf)!"
        gzip -c $high_pos.temp > $raw_high
    fi
   
    gunzip -c $raw_high > raw_high.temp
    cat $high_pos.temp raw_high.temp > raw_high_updated.temp
    gzip -c raw_high_updated.temp > $raw_high
  
    gzip -c $high_pos.temp > $high_pos

    # --- low-confidence --- #
    if ! [ -a $raw_low ]
    then
        echo "Creating master list (low conf)!"
        gzip -c $low_pos.temp > $raw_low
    fi

    gunzip -c $raw_low > raw_low.temp
    cat $low_pos.temp raw_low.temp > raw_low_updated.temp
    gzip -c raw_low_updated.temp > $raw_low

    gzip -c $low_pos.temp > $low_pos

done

gunzip -c $raw_high | sort-bed - | bedops --merge - | awk '{ print $1, $2, $3, "FP_NA_"NR }' > merged_high.temp
gunzip -c $raw_low | sort-bed - | bedops --merge - | awk '{ print $1, $2, $3, "FP_lNA_"NR }' > merged_low.temp

# --- now check overlaps with the time-points, to find out what contributed and when... --- #
for datatime in 0min 2min 5min 12.5min 25min 50min
do
    echo "Overlapping with" $datatime
    high_pos=$results/$datatype\_$datatime.pos.$THRE_high.bedGraph.gz
    low_pos=$results/$datatype\_$datatime.pos.$THRE_low.bedGraph.gz

    bedmap --echo --echo-map --multidelim '|' merged_high.temp $high_pos.temp > merged_high2.temp
    mv merged_high2.temp merged_high.temp
    bedmap --echo --echo-map --multidelim '|' merged_low.temp $low_pos.temp >  merged_low2.temp
    mv merged_low2.temp merged_low.temp
done

# Call on the other script to tidy this up!
echo "Master files created, processing!"
python process_master.py merged_high.temp $proc_high $datatype
python process_master.py merged_low.temp $proc_low $datatype

echo "Looking at overlap with genomic regions and ChIPseq data!"
awk 'BEGIN{OFS="\t"}{ if ($6=="+") { print $1, $2, $2+1, $4, $6 } else { print $1, $3-1, $3, $4, $6 } }' $gene_list > gene_starts.temp

# (stripping the header line)
gunzip -c $proc_high | sed '/inflation/d' > pred_high.temp
gunzip -c $proc_low | sed '/inflation/d' > pred_low.temp

# Find out where the hits are!
echo "Dividing hits into regions."
bedmap --fraction-ref 0.5 --indicator pred_high.temp $gene_start > gene_start_high.temp
bedmap --fraction-ref 0.5 --indicator pred_low.temp $gene_start > gene_start_low.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred_high.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb_high.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred_low.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb_low.temp
bedmap --fraction-ref 0.5 --indicator pred_high.temp $non_gene > non_gene_high.temp
bedmap --fraction-ref 0.5 --indicator pred_low.temp $non_gene > non_gene_low.temp

# Get the distance to the nearest gene_start... (for all hits!) (modify the master)
echo "Getting distance to closest gene start."
bedtools closest -d -a pred_high.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-2), $NF, $(NF-1) }' | python post_bedtools_join_closest.py > dist_high.temp
bedtools closest -d -a pred_low.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-2), $NF, $(NF-1) }' | python post_bedtools_join_closest.py > dist_low.temp
# (creates dist_high.temp and dist_low.temp)

# Get histone counts! (totally unnormalised here...)
echo "Getting histone counts."
gunzip -c $H3K4me1 | bedmap --delim '\t' --range 500 --sum pred_high.temp - > H3K4me1_high.temp
gunzip -c $H3K4me3 | bedmap --delim '\t' --range 500 --sum pred_high.temp - > H3K4me3_high.temp
gunzip -c $H3K27ac | bedmap --delim '\t' --range 500 --sum pred_high.temp - > H3K27ac_high.temp
gunzip -c $H3K4me1 | bedmap --delim '\t' --range 500 --sum pred_low.temp - > H3K4me1_low.temp
gunzip -c $H3K4me3 | bedmap --delim '\t' --range 500 --sum pred_low.temp - > H3K4me3_low.temp
gunzip -c $H3K27ac | bedmap --delim '\t' --range 500 --sum pred_low.temp - > H3K27ac_low.temp

echo "Combining!"
# create a headerfile
gunzip -c $proc_high | head -1 | awk '{ print $0, "gene_start","gene_body","non_gene","closest_gene","distance","strand","H3K4me1","H3K4me3","H3K27ac" }' > headerfile.temp
# stick it all together
paste pred_high.temp gene_start_high.temp gb_high.temp non_gene_high.temp dist_high.temp H3K4me1_high.temp H3K4me3_high.temp H3K27ac_high.temp > comb_high.temp
cat headerfile.temp comb_high.temp | gzip -c > $pre_QC_confident

paste pred_low.temp gene_start_low.temp gb_low.temp non_gene_low.temp dist_low.temp H3K4me1_low.temp H3K4me3_low.temp H3K27ac_low.temp > comb_low.temp
cat headerfile.temp comb_low.temp | gzip -cv > $maybe

echo "Doing quality control!"
R --slave --file=QC.r --args $pre_QC_confident m.temp
gzip -c m.temp > $confident

echo "Dividing results into regions, for visualisation." # (this used to be 'extract_regions'
awk 'BEGIN{OFS="\t"}{ if ($13==1) print $1,$2,$3,$4,$16,$18 }' m.temp | sed '1d' | sed '/+;-/d' > $results/metaplot/dREG_gs.bed
awk 'BEGIN{OFS="\t"}{ if (($14!=0)&&($17>2000)) print $1,$2,$3,$4,$16,$18 }' m.temp | sed '1d'| sed '/+;-/d'  > $results/metaplot/dREG_gb.bed
# non-genes don't have strands!
awk 'BEGIN{OFS="\t"}{ if ($15==1) print $1,$2,$3,$4,$16 }' m.temp | sed '1d' | sed '/+;-/d' > $results/metaplot/dREG_ng.bed

# -- now for low... (will remove this after sanity checks, I think) ...
awk 'BEGIN{OFS="\t"}{ if ($13==1) print $1,$2,$3,$4,$16,$18 }' comb_low.temp| sed '/+;-/d' > $results/metaplot/dREG_maybe_gs.bed
awk 'BEGIN{OFS="\t"}{ if (($14!=0)&&($17>2000)) print $1,$2,$3,$4,$16,$18 }' comb_low.temp | sed '/+;-/d' > $results/metaplot/dREG_maybe_gb.bed
awk 'BEGIN{OFS="\t"}{ if ($15==1) print $1,$2,$3,$4,$16 }' comb_low.temp  | sed '/+;-/d' > $results/metaplot/dREG_maybe_ng.bed

# --- Create a log file of what just happened --- #
echo `date` > $results/logfile.txt
echo "Confident threshold:" $THRE_high >> $results/logfile.txt
echo "Maybe threshold:" $THRE_low >> $results/logfile.txt
echo "Datatype:" $datatype >> $results/logfile.txt
echo "Suffix:" $suffix >> $results/logfile.txt
echo "# High confidence regions:" `wc -l comb_high.temp | awk '{ print $1 }'` >> $results/logfile.txt
echo "# Low confidence regions:" `wc -l comb_low.temp | awk '{ print $1 }'` >> $results/logfile.txt

# --- Clean up ---- #
echo "Tidying up!"
#rm -v *.temp
#rm -v $raw_high
#rm -v $raw_low
#rm -v $proc_high
#rm -v $proc_low
