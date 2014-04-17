#!/bin/bash
# author:       Stephanie Hyland (sh985@cornell.edu)
# date:         February 2014
# description:  Takes dREG output, thresholds into positive regions.
#               Adds new regions to master list.

# Threshold
THRE=0.5

# Inputs!
datatype=$1     # eg FP, TRP, control
#datatime=$2     # eg 2min, 5min, etc... or DMSO

# Paths
data=/Users/stephanie/ll/data/
results=/Users/stephanie/ll/results
idlist=$results/$datatype\_dREG_regions.bed.gz

for datatime in 0min 2min 5min 12.5min 25min 50min
do
    echo "Encorporating" $datatime
    # Set up vars
    folder=$results/$datatype/$datatime
    preds=$folder/$datatype\_$datatime.predictions.bedGraph.gz
    pos=${preds/predictions/positive}

    # What is the score distribution? (can plot this later as necessary)
    echo "Getting scores."
    gunzip -c $preds | awk '{ print $4 }' > ${preds/.predictions.bedGraph.gz/.scores.txt}
    # Expands positive regions into window +-50, merges these.
    echo "Getting regions."
    gunzip -c $preds | awk 'BEGIN{OFS="\t"} ($4 > '"$THRE"') { print $1,$2-50,$3+51,$4 }' | sort-bed - | bedops --merge - | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, "'$datatype'_" "'${datatime/min/}'_" NR }' > $pos.temp

    # Create the master list, if necessary
    if ! [ -a $idlist ]
    then
        echo "Creating master list!"
        gzip -c $pos.temp > $idlist
    fi

    # Checks if this data already exists in the master list
    if gunzip -c $idlist | grep -q $datatype\_${datatime/min/}_
    then
        echo "Master list already includes these regions - skipping."
    else
        # Adds regions to existing list (assume same if ANY overlap between sets)
        echo "Updating master list."
        gunzip -c $idlist | bedmap --delim '|' --multidelim '|' --echo --indicator $pos.temp - > new_or_old.temp
        awk 'BEGIN{FS="|"}{ if ($2==0) print $1 }' new_or_old.temp > new.temp
        awk 'BEGIN{FS="|"}{ if ($2==1) print $1 }' new_or_old.temp > old.temp
        gunzip -c $idlist | bedmap --delim '|' --multidelim '|' --echo --echo-map - old.temp  > idlist.temp
        cat idlist.temp new.temp | sort-bed - > idlist_updated.temp
        gzip -c idlist_updated.temp > $idlist
    fi

    # Tidy up
    gzip -c $pos.temp > $pos
    echo "Deleting files."
    rm -v *.temp
done

# Call on the other script to tidy this up!
echo "Master file created, processing!"
python process_master.py $datatype

echo "Looking at overlap with genomic regions and ChIPseq data!"
./get_regions.sh ${idlist/.bed.gz/_uniq.bed.gz}
