#!/bin/bash
# Convert bed files to bigwig - separating strands, of course!

mm9=/Users/stephanie/ll/data/mm9.chromSizes

for f in *.bed
do
    echo $f
    echo "Separating strands."
    python ~/Projects/feabhsa-here/sep_strand.py $f

    plus=$f.plus
    minus=$f.minus

    echo "Converting to bedGraph!"
    /Users/stephanie/Tools/UCSC/bedItemOverlapCount -chromSize=$mm9 mm9 $plus > ${f/.bed/.plus.bedGraph}
    /Users/stephanie/Tools/UCSC/bedItemOverlapCount -chromSize=$mm9 mm9 $minus > ${f/.bed/.minus.bedGraph}

    echo "Converting to bigWig!"
    /Users/stephanie/Tools/UCSC/bedGraphToBigWig ${f/.bed/.plus.bedGraph} $mm9 ${f/.bed/.plus.bw}
    /Users/stephanie/Tools/UCSC/bedGraphToBigWig ${f/.bed/.minus.bedGraph} $mm9 ${f/.bed/.minus.bw}

    echo "Compressing!"
    gzip -v *.bedGraph

    echo "Tidying!"
    rm -v $f.plus
    rm -v $f.minus
done
