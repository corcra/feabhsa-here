#!/bin/bash

#suffix=full_SVM      # eg full_SVM, 40m_SVM...
suffix=40m_SVM      # eg full_SVM, 40m_SVM...
#version=C           # eg B, C, V6.5, etc...
#version=B           # eg B, C, V6.5, etc...
version=R           # eg B, C, V6.5, etc...
#version=V6.5
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
echo "Ran HMM with options: -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt" >> $logfile
./code/hmm2 -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt

# --- do some QC! --- #
echo "Starting QC..."
echo `date` > lost_in_QC.$version.$time.txt
# add back the gene strand info ... resulting columns are NAME, GID, LEN, ROUNDS, TRANSITION, DENSITY1, DENSITY2, STRAND, CHR, START, END
sort -k2,2 $result > res_sort.$version.$time.temp
sort -k2,2 $genelist | awk '{ print $2, $5, $4, $6, $7 }' | join -1 2 -2 1 res_sort.$version.$time.temp - > res.$version.$time.temp

# rounds = 200, kick it out
awk '{ if ($4!=200) print $0 }' res.$version.$time.temp > res2.$version.$time.temp
echo $[`wc -l res.$version.$time.temp | awk '{ print $1 }'` - `wc -l res2.$version.$time.temp | awk '{print $1}'`] "genes had >200 rounds - removed." >> $logfile
awk '{ if ($4==200) print $0,"rounds" }' res.$version.$time.temp >> lost_in_QC.$version.$time.txt

# transition > len, kick it out
awk '{ if ($5<$3) print $0 }' res2.$version.$time.temp > res3.$version.$time.temp
echo $[`wc -l res2.$version.$time.temp | awk '{ print $1 }'` - `wc -l res3.$version.$time.temp | awk '{print $1}'`] "genes had transition > len - removed." >> $logfile
awk '{ if ($5>=$3) print $0, "transition>len" }' res2.$version.$time.temp >> lost_in_QC.$version.$time.txt

# transition == 2*binsize, kick it out
awk '{ if ($5!=2*'$binsize') print $0 }' res3.$version.$time.temp > res4.$version.$time.temp
echo $[`wc -l res3.$version.$time.temp | awk '{ print $1 }'` - `wc -l res4.$version.$time.temp | awk '{print $1}'`] "genes had transition == 2*binsize - removed." >> $logfile
awk '{ if ($5==2*'$binsize') print $0,"transition=2binsize" }' res3.$version.$time.temp >> lost_in_QC.$version.$time.txt

# density1 > density2, kick it out
awk '{ if ($6<0.5*$7) print $0 }' res4.$version.$time.temp > res5.$version.$time.temp
echo $[`wc -l res4.$version.$time.temp | awk '{ print $1 }'` - `wc -l res5.$version.$time.temp | awk '{print $1}'`] "genes had density1 > 0.5*density2 - removed." >> $logfile
awk '{ if ($6>=0.5*$7) print $0,"density" }' res4.$version.$time.temp >> lost_in_QC.$version.$time.txt

# transition involves scientific notation, kick it out
grep -v "e" res5.$version.$time.temp > res6.$version.$time.temp
echo $[`wc -l res5.$version.$time.temp | awk '{ print $1 }'` - `wc -l res6.$version.$time.temp | awk '{print $1}'`] "genes had scientific notation in their transition - removed." >> $logfile
grep "e" res5.$version.$time.temp | awk '{ print $0, "sci-notation" }' >> lost_in_QC.$version.$time.txt

# transition overlap with a dREG hit? kick it out
# make beddy, e.g. CHR, TRANSITION_START, TRANSITION_END, STRAND, GID, NAME, LEN, ROUNDS, TRANSITION, DENSITY1, DENSITY2, GENE_START, GENE_END
awk '{ if ($8=="+") {{ printf "%s %i %i ", $9, $10+$5, $10+$5+1} { print $8, $1, $2, $3, $4, $5, $6, $7, $10, $11 }} else  {{ printf "%s %i %i ", $9, $11-$3-1, $11-$3} { print $8, $1, $2, $3, $4, $5, $6, $7, $10, $11 }}}' res6.$version.$time.temp | sort-bed - > res.bed.$version.$time.temp
gunzip -c $dREG_list | sed '1d' | bedmap --range $binsize --echo --indicator res.bed.$version.$time.temp - > res2.bed.$version.$time.temp
grep '|0' res2.bed.$version.$time.temp | awk 'BEGIN{FS="|"}{print $1}' > res3.bed.$version.$time.temp
echo $[`wc -l res.bed.$version.$time.temp | awk '{ print $1 }'` - `wc -l res3.bed.$version.$time.temp | awk '{print $1}'`] "genes had a dREG hit near their transition - removed." >> $logfile
grep '|1' res2.bed.$version.$time.temp | awk 'BEGIN{FS="|"}{print $1 }' | awk '{print $5, $6, $7, $8, $9, $10, $11, $4, $1, $12, $13, "dREG" }' >> lost_in_QC.$version.$time.txt

echo "Overall," $[`wc -l res.$version.$time.temp | awk '{ print $1 }'` - `wc -l res3.bed.$version.$time.temp | awk '{print $1}'`] "genes removed for QC." >> $logfile

# convert it back into something nice (by nice I really just mean 'something along the lines of the output of the HMM'...)
echo "gid" "name" "len" "rounds" "transition" "density1" "density2" "strand" "chr" "gene_start" "gene_end" > final.$version.$time.temp
awk 'BEGIN{OFS="\t"}{ print $6, $5, $7, $8, $9, $10, $11, $4, $1, $12, $13 }' res3.bed.$version.$time.temp >> final.$version.$time.temp

echo "After QC, there are" `wc -l final.$version.$time.temp| awk '{ print $1 }'` "genes remaining." >> $logfile
mv -v final.$version.$time.temp $result

# --- Tidy! --- #
rm *$version.$time.temp

# --- Format the deleted list --- #
awk '{ print $9, $10, $11, $2, $1, $8, $3, $4, $5, $6, $7, $12 }' lost_in_QC.$version.$time.txt | sed '1d' | sort-bed - > lost_in_QC.$version.$time.bed
rm -v lost_in_QC.$version.$time.txt
