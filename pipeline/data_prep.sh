#!/bin/bash
# author:       Stephanie Hyland (sh985@cornell.edu)
# date:         February-April 2014
# description:  This is something of a 'master' pipeline. It does the following:
#               1. Combine dREG predictions (above thresh) from all timepoints (into $rawlist).
#               2. Format this nicely (via process_master.py)
#               3. Classify all points by genomic region, get overlap with histone marks.

# Threshold
THRE=0.8

# Inputs!
datatype=FP     # eg FP, TRP... not sure why I'm still asking this
suffix=full_SVM         # eg full_SVM, 40m_SVM... other?

# Paths
data=/Users/stephanie/ll/data
folder=$data/dREG/$suffix
results=/Users/stephanie/ll/results/$suffix

# Main files
rawlist=$results/$datatype\_dREG_regions.bed.gz
regionlist=$results/$datatype\_dREG_regions_uniq.bed.gz
final=$results/$datatype\_dREG_regions_marked.bed.gz

# Genes (for later)
gene_list=$data/genes/gene_list.bed
gene_start=$data/genes/gene_start.bed
gene_body=$data/genes/gene_body.bed
non_gene=$data/genes/non_genes.bed

# Markers (also for later)
H3K4me1=$data/ChIPseq/H3K4me1_Meiss.sorted.bed.gz
H3K4me3=$data/ChIPseq/H3K4me3_Marson.sorted.bed.gz
H3K27ac=$data/ChIPseq/H3K27Ac.sorted.bed.gz

for datatime in 0min 2min 5min 12.5min 25min 50min
do
    echo "Encorporating" $datatime
    # Set up vars
    preds=$folder/$datatype\_$datatime.predictions.bedGraph.gz
    pos=$results/$datatype\_$datatime.positive.bedGraph.gz

    # What is the score distribution? (can plot this later as necessary)
    echo "Getting scores."
    gunzip -c $preds > $preds.temp
    bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $gene_start | gzip -c > ${preds/.predictions/.gs.predictions}
    bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $gene_body | gzip -c > ${preds/.predictions/.gb.predictions}
    bedmap --fraction-ref 0.5 --echo --skip-unmapped $preds.temp $non_gene | gzip -c > ${preds/.predictions/.ng.predictions}

    # Expands positive regions into window +-50, merges these.
    echo "Getting regions."
    awk 'BEGIN{OFS="\t"} ($4 > '"$THRE"') { print $1,$2-50,$3+51,$4 }' $preds.temp | sort-bed - | bedops --merge - | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, "'$datatype'_" "'${datatime/min/}'_" NR }' > $pos.temp

    # Create the master list, if necessary
    if ! [ -a $rawlist ]
    then
        echo "Creating master list!"
        gzip -c $pos.temp > $rawlist
    fi

    # Checks if this data already exists in the master list
    if gunzip -c $rawlist | grep -q $datatype\_${datatime/min/}_
    then
        echo "Master list already includes these regions - skipping."
    else
        # Adds regions to existing list (assume same if ANY overlap between sets)
        echo "Updating master list."
        gunzip -c $rawlist | bedmap --delim '|' --multidelim '|' --echo --indicator $pos.temp - > new_or_old.temp
        awk 'BEGIN{FS="|"}{ if ($2==0) print $1 }' new_or_old.temp > new.temp
        awk 'BEGIN{FS="|"}{ if ($2==1) print $1 }' new_or_old.temp > old.temp
        gunzip -c $rawlist | bedmap --delim '|' --multidelim '|' --echo --echo-map - old.temp  > rawlist.temp
        cat rawlist.temp new.temp | sort-bed - > rawlist_updated.temp
        gzip -c rawlist_updated.temp > $rawlist
    fi

    # Tidy up
    gzip -c $pos.temp > $pos
    echo "Deleting files."
    rm -v *.temp
    rm -v $pos.temp
    rm -v $preds.temp
done

# Call on the other script to tidy this up!
echo "Master file created, processing!"
python process_master.py $results $datatype

echo "Looking at overlap with genomic regions and ChIPseq data!"

awk 'BEGIN{OFS="\t"}{ if ($6=="+") { print $1, $2, $2+1, $4, $6 } else { print $1, $3-1, $3, $4, $6 } }' $gene_list > gene_starts.temp

# You know what, just give me the whole path.
gunzip -c $regionlist | sed '/inflation/d' > pred.temp
# (stripping the header line)

# Find out where the hits are!
echo "Dividing hits into regions."
bedmap --fraction-ref 0.5 --indicator pred.temp $gene_start > gene_start.temp
bedmap --fraction-ref 0.5 --delim ';' --multidelim ';' --echo-map-id pred.temp $gene_body | awk '{ if (NF==0) print 0; else print $0 }' > gb.temp
bedmap --fraction-ref 0.5 --indicator pred.temp $non_gene > non_gene.temp

# Get the distance to the nearest gene_start... (for all hits!) (modify the master)
echo "Getting distance to closest gene start."
bedtools closest -d -a pred.temp -b gene_starts.temp | awk 'BEGIN{OFS="\t"}{ print $4, $(NF-2), $NF, $(NF-1) }' > dist.pre.temp
python post_bedtools_join_closest.py
# (creates dist.temp)

# Get histone counts! (totally unnormalised here...)
echo "Getting histone counts."
gunzip -c $H3K4me1 | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K4me1.temp
gunzip -c $H3K4me3 | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K4me3.temp
gunzip -c $H3K27ac | bedmap --delim '\t' --range 500 --sum pred.temp - > H3K27ac.temp

echo "Combining!"
paste pred.temp gene_start.temp gb.temp non_gene.temp dist.temp H3K4me1.temp H3K4me3.temp H3K27ac.temp > comb.temp
gunzip -c $regionlist | head -1 | awk '{ print $0, "gene_start","gene_body","non_gene","closest_gene","distance","strand","H3K4me1","H3K4me3","H3K27ac" }' > headerfile.temp
cat headerfile.temp comb.temp | gzip -c > $final

echo "Tidying up!"
rm -v *.temp
