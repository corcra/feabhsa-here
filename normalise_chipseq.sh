#!/bin/bash
# Normalise histone mods (H3) by windows of total H3, with pseudocounts

windowsize=250
pseudocount=1

# Chipseq files (bedGraphs!)
me1=/Users/stephanie/ll/data/ChIPseq/bedgraphs/H3K4me1_Meiss.sorted.bedGraph.gz
me3=/Users/stephanie/ll/data/ChIPseq/bedgraphs/H3K4me3_Marson.sorted.bedGraph.gz
ac=/Users/stephanie/ll/data/ChIPseq/bedgraphs/H3K27Ac.sorted.bedGraph.gz

# note: this is a bed file...
H3_path=/Users/stephanie/ll/data/ChIPseq/H3_align.sorted.bed

# genome info (for converting to bw)
mm9=/Users/stephanie/ll/data/mm9.chromSizes

echo "Creating fake windows."
python fake_windows.py $windowsize

echo "Getting counts in windows."
cat $H3_path | bedmap --delim '\t' --echo --count windows_$windowsize.bed - > windows_$windowsize\_counts.bed

echo "Applying normalisation."
gunzip -c $me1 | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me1/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $me3 | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me3/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $ac | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${ac/.bedGraph.gz/.window$windowsize\norm.bedGraph} 

echo "Converting to bw!"
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me1/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me1/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me3/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me3/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${ac/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${ac/.bedGraph.gz/.window$windowsize\norm.bw}

mv -v /Users/stephanie/ll/data/ChIPseq/bedgraphs/*.bw /Users/stephanie/ll/data/ChIPseq/bw/
