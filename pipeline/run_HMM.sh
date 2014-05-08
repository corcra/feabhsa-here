#!/bin/bash

#suffix=full_SVM      # eg full_SVM, 40m_SVM...
suffix=40m_SVM      # eg full_SVM, 40m_SVM...
version=B           # eg B, C, V6.5, etc...
THRE_high=0.8

hmm_dir=/Users/stephanie/ll/from_hojoong
cd $hmm_dir

bg_dir=/Users/stephanie/ll/data/FP/bedgraphs
time=$1

dREG_list=/Users/stephanie/ll/results/$suffix/dREG_regions_confident_$THRE_high.bed.gz
genes_preformatted=/Users/stephanie/ll/results/$suffix/active_genes.bed

if [ $time == "5" ]
then
    thre=30000
    binsize=500
fi

if [ $time == "12.5" ] || [ $time == "1230" ]
then
    thre=30000
    binsize=1000
fi

if [ $time == "25" ]
then
    thre=60000
    binsize=2000
fi

if [ $time == "50" ]
then
    thre=150000
    binsize=5000
fi

genelist=$hmm_dir/list/genes_over$thre.txt
logfile=$hmm_dir/result/$version\_$time\min_bin$binsize.$suffix.logfile.txt

echo `date` > $logfile
echo "Suffix is" $suffix >> $logfile
echo "Using threshold of" $thre "and binsize of" $binsize "because time is" $time"." >> $logfile

bg_plus=$bg_dir/$version_$1\minFP_Plus.bedGraph
bg_minus=$bg_dir/$version_$1\minFP_Minus.bedGraph

bi_plus=$hmm_dir/bi/$version\_$time\min.p.bi
bi_minus=$hmm_dir/bi/$version\_$time\min.m.bi
bi0_plus=$hmm_dir/bi/$version\_0min.p.bi
bi0_minus=$hmm_dir/bi/$version\_0min.m.bi

result=$hmm_dir/result/$version\_$time\min_bin$binsize.$suffix.txt

# --- prepare the genelist --- #
if ! [ -a $genelist ]
    then
        echo "Processing genelist!"
        awk 'BEGIN{OFS="\t"}{ if (($3-$2)>'$thre') print NR,$4,$4,$1,$6,$2,$3 }' $genes_preformatted > $genelist
        echo "There are" `wc -l $genelist | awk '{print $1 }'` "genes remaining!"
fi

echo `wc -l $genelist | awk '{print $1 }'` "starting genes to analyse." >>$logfile

#  --- prepare the bedgraphs --- #
if ! [ -e $bi_plus ]
    then
        if [ -e $bg_plus.gz ] && ! [ -e $bg_plus ]
            then
                echo "Unzipping bedgraph ("$bg_plus")"
                gunzip $bg_plus.gz
        fi
        echo "Creating index for" $bg_plus
        ./code/prereadbg -i $bg_plus -o $bi_plus -c list/mm9chr.txt
fi

if ! [ -e $bi_minus ]
    then
        if [ -e $bg_minus.gz ] && ! [ -e $bg_minus ]
            then
                echo "Unzipping bedgraph ("$bg_minus")"
                gunzip $bg_minus.gz
        fi
        echo "Creating index for" $bg_minus
        ./code/prereadbg -i $bg_minus -o $bi_minus -c list/mm9chr.txt
fi

# --- run the HMM --- #
echo "Running HMM!"
echo "Options: -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt"
./code/hmm2 -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt

# --- check for gene-body enhancers --- #
# Have a file of hits specifically in the gene body...
# sort the genelist
sort -k4,4 $genes_preformatted > genes.temp
# sort the result-list, then perform a join on the genename field... the result should be the addition of two columns - first is the number of dREG hits in the gene start region, second is number of hits in the gene body ... we can use this for later analysis!
sort -k2,2 $result | join -1 2 -2 4 - genes.temp  > res.temp

# --- do some QC! --- #
# rounds = 200, kick it out
awk '{ if ($4!=200) print $0 }' res.temp > res2.temp
echo $[`wc -l res.temp | awk '{ print $1 }'` - `wc -l res2.temp | awk '{print $1}'`] "genes had >200 rounds - removed." >> $logfile
# transition > len, kick it out
awk '{ if ($5<$3) print $0 }' res2.temp > res3.temp
echo $[`wc -l res2.temp | awk '{ print $1 }'` - `wc -l res3.temp | awk '{print $1}'`] "genes had transition > len - removed." >> $logfile
# transition == 2*binsize, kick it out
awk '{ if ($5!=2*'$binsize') print $0 }' res3.temp > res4.temp
echo $[`wc -l res3.temp | awk '{ print $1 }'` - `wc -l res4.temp | awk '{print $1}'`] "genes had transition == 2*binsize - removed." >> $logfile
# density1 > density2, kick it out
awk '{ if ($6<0.5*$7) print $0 }' res4.temp > res5.temp
echo $[`wc -l res4.temp | awk '{ print $1 }'` - `wc -l res5.temp | awk '{print $1}'`] "genes had density1 > 0.5*density2 - removed." >> $logfile
# transition overlap with a dREG hit? kick it out
awk '{ if (12=="+") { print $8, $9+3,$9+3+1,$1,$12,$2,$3,$4,$5,$6,$7,$9,$10,$13,$14} else { print $8, $10-$3-1, $10-$3, $1,$12,$2,$3,$4,$5,$6,$7,$9,$10,$13,$14} }' res5.temp | sort-bed - > res.bed.temp
gunzip -c $dREG_list | sed '1d' | bedmap --range $binsize --echo --indicator res.bed.temp - | grep '|0' | awk 'BEGIN{FS="|"}{print $1}' > res2.bed.temp
echo $[`wc -l res.bed.temp | awk '{ print $1 }'` - `wc -l res2.bed.temp | awk '{print $1}'`] "genes had a dREG hit near their transition - removed." >> $logfile

echo "Overall," $[`wc -l res.temp | awk '{ print $1 }'` - `wc -l res2.bed.temp | awk '{print $1}'`] "genes removed for QC." >> $logfile
# convert it back into something nice (by nice I really just mean 'something along the lines of the output of the HMM'...)
echo "gid" "name" "len" "rounds" "transition" "density1" "density2" "chr" "gene_start" "gene_end" "n_gs" "n_gb" > final.temp
awk 'BEGIN{OFS="\t"}{ print $6, $4, $7, $8, $9, $10, $11, $1, $12, $13, $14, $15 }' res2.bed.temp >> final.temp

echo "After QC, there are" `wc -l final.temp| awk '{ print $1 }'` "genes remaining." >> $logfile
mv final.temp $result

# --- Tidy! --- #
rm *.temp
