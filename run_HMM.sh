#!/bin/bash

suffix=40m_SVM      # eg full_SVM, 40m_SVM...

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

genes_preformatted=/Users/stephanie/ll/results/$suffix/genes_gs_hit.bed
genelist=$hmm_dir/list/genes_over$thre.txt

result=$hmm_dir/result/V6.5_$time\min_bin$binsize.txt

# --- prepare the genelist --- #
if ! [ -a $genelist ]
    then
        echo "Processing genelist!"
        awk 'BEGIN{OFS="\t"}{ if ((($3-$2)>'$thre')&&($7>0)) print NR,$4,$4,$1,$6,$2,$3 }' $genes_preformatted > $genelist
        echo "There are" `wc -l $genelist | awk '{print $1 }'` "genes remaining!"
fi


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
./code/hmm2 -o $result -p $bi_plus -m $bi_minus -p0 $bi0_plus -m0 $bi0_minus -b $binsize -g $genelist -c list/mm9chr.txt

# --- check for gene-body enhancers --- #
# Have a file of hits specifically in the gene body...
echo `sed '1q;d' $result` "chr" "start" "end" "n_gs" "n_gb" > header.temp
# sort the genelist
sort -k4,4 $genes_preformatted > genes.temp
# sort the result-list, then perform a join on the genename field... the result should be the addition of two columns - first is the number of dREG hits in the gene start region, second is number of hits in the gene body ... we can use this for later analysis!
sort -k2,2 $result | join -1 2 -2 4 - genes.temp | awk 'BEGIN{OFS="\t"}{ print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }' > $result.temp
cat header.temp $result.temp > $result.final
mv $result.final $result

# --- Tidy! --- #
rm -v *.temp
