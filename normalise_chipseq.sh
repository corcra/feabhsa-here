#!/bin/bash
# Normalise histone mods (H3) by windows of total H3, with pseudocounts

windowsize=1000
pseudocount=1

# Chipseq files (bedGraphs!)
me1_plus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K4me1_Meiss.sorted.plus.bedGraph.gz
me1_minus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K4me1_Meiss.sorted.minus.bedGraph.gz
me3_plus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K4me3_Marson.sorted.plus.bedGraph.gz
me3_minus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K4me3_Marson.sorted.minus.bedGraph.gz
ac_plus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K27Ac.sorted.plus.bedGraph.gz
ac_minus=/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep/H3K27Ac.sorted.minus.bedGraph.gz

# note: this is a bed file...
H3_path=/Users/stephanie/ll/data/ChIPseq/H3_align.sorted.bed

# genome info (for converting to bw)
mm9=/Users/stephanie/ll/data/mm9.chromSizes

echo "Creating fake windows."
python fake_windows.py $windowsize

echo "Getting counts in windows."
sed -n '/+/p' $H3_path | bedmap --delim '\t' --echo --count windows_$windowsize.bed - > windows_$windowsize\_counts.plus.bed
sed -n '/-/p' $H3_path | bedmap --delim '\t' --echo --count windows_$windowsize.bed - > windows_$windowsize\_counts.minus.bed

echo "Applying normalisation."
gunzip -c $me1_plus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.plus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me1_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $me1_minus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.minus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me1_minus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $me3_plus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.plus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me3_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $me3_minus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.minus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${me3_minus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $ac_plus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.plus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${ac_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 
gunzip -c $ac_minus | bedmap --delim '\t' --echo --echo-map-score - windows_$windowsize\_counts.minus.bed | awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4/($5+'$pseudocount') }'  > ${ac_minus/.bedGraph.gz/.window$windowsize\norm.bedGraph} 

echo "Converting to bw!"
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me1_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me1_plus/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me1_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me1_minus/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me3_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me3_plus/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${me3_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${me3_minus/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${ac_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${ac_plus/.bedGraph.gz/.window$windowsize\norm.bw}
/Users/stephanie/Tools/UCSC/bedGraphToBigWig ${ac_plus/.bedGraph.gz/.window$windowsize\norm.bedGraph} $mm9 ${ac_minus/.bedGraph.gz/.window$windowsize\norm.bw}
