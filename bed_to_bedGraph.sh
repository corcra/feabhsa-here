#!/bin/bash
# Converts a bed to a bedGraph, 
# THEN NORMALISES BY READ COUNT (not really!)

#mm9=/home/macproadmin/users/Stephanie/data/mm9.chromSizes
mm9=/Users/stephanie/ll/data/mm9.chromSizes

for f in *.bed
do
    echo $f
    echo "Separating strands."
    python ~/Projects/feabhsa-here/sep_strand.py $f

    echo "Getting read depth."
    plus=$f.plus
    n_plus=`wc -l $plus | cut -d ' ' -f 1`
    minus=$f.minus
    n_minus=`wc -l $minus | cut -d ' ' -f 1`
    
    echo $plus
    echo $n_plus
    echo $minus
    echo $n_minus
    read -p "Press [Enter] key to continue..."

    echo "Converting to bedGraph!"
    /Users/stephanie/Tools/UCSC/bedItemOverlapCount -chromSize=$mm9 mm9 $plus > plus.bedGraph
    /Users/stephanie/Tools/UCSC/bedItemOverlapCount -chromSize=$mm9 mm9 $minus > minus.bedGraph

#    echo "Normalising by total reads!"
#    awk '{ print $1, $2, $3, $4/$n_plus }' plus.bedGraph > plus.norm
#    awk '{ print $1, $2, $3, $4/$n_minus }' minus.bedGraph > minus.norm

    echo "Combining strands and sorting!"
    cat plus.norm minus.norm | sort-bed - > ${f/.bed/.bedGraph}

    echo "Tidying!"
    rm -v plus.bedGraph
    rm -v minus.bedGraph
    rm -v plus.norm
    rm -v minus.norm

done
