#!/bin/bash

hmm_dir=/Users/stephanie/ll/from_hojoong
cd $hmm_dir

bg_dir=/Users/stephanie/ll/data/FP/bedgraphs
time=$1
binsize=$2
#thre=$((2*$binsize))
thre=25000
#if [ $thre -lt $((2*$binsize)) ]
#then
#    echo "Genes must be at least twice as long as the binsize! Modifying. New threshold:" $thre
#fi
#thre=1000

bg_plus=$bg_dir/V6.5_$1\minFP_Plus.bedGraph
bg_minus=$bg_dir/V6.5_$1\minFP_Minus.bedGraph

bi_plus=$hmm_dir/bi/V6.5_$time\min.p.bi
bi_minus=$hmm_dir/bi/V6.5_$time\min.m.bi
bi0_plus=$hmm_dir/bi/V6.5_0min.p.bi
bi0_minus=$hmm_dir/bi/V6.5_0min.m.bi

genes_preformatted=/Users/stephanie/ll/data/genes/analyse_FP/gene_list_with_overlaps.bed
genelist=$hmm_dir/list/genes_over$thre.txt

result=$hmm_dir/result/V6.5_$time\min_bin$binsize.txt

# prepare the genelist
if ! [ -a $genelist ]
    then
        echo "Processing genelist!"
        awk 'BEGIN{OFS="\t"}{ if ((($3-$2)>'$thre')&&($7>0)) print NR,$4,$4,$1,$6,$2,$3 }' $genes_preformatted > $genelist
fi

echo "There are" `wc -l $genelist | awk '{print $1 }'` "genes remaining!"

# prepare the bedgraphs
if ! [ -a $bi_plus ]
    then
        if ! [ -a $bg_plus.gz ]
            then
                echo "Unzipping bedgraph ("$bg_plus")"
                gunzip $bg_plus.gz
        fi
        echo "Creating index for" $bg_plus
        ./code/prereadbg -i $bg_plus -o $bi_plus -c list/mm9chr.txt
fi

if ! [ -a $bi_minus ]
    then
        if ! [ -a $bg_minus.gz ]
            then
                echo "Unzipping bedgraph ("$bg_minus")"
                gunzip $bg_minus.gz
        fi
        echo "Creating index for" $bg_minus
        ./code/prereadbg -i $bg_minus -o $bi_minus -c list/mm9chr.txt
fi

# run the HMM
./code/hmm2 -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt
