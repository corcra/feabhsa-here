#!/usr/bin/python
# just make a fake file of windows (over entire mouse genome) (we can use this to make normalisation factor in each window!)

import sys

chro_data='/Users/stephanie/ll/data/mm9.chromSizes'
window_size=int(sys.argv[1])

chro_sizes_dat=open(chro_data).readlines()
chro_sizes_list=[(chro.split()[0],int(chro.split()[1])) for chro in chro_sizes_dat]
chro_sizes={chro:size for (chro,size) in chro_sizes_list}

fakefile=open('windows_'+str(window_size)+'.bed','w')

# this is the order chromosomes appear in sorted bed files...
chro_list=['chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr8','chr9','chrM','chrX','chrY']

for chro in chro_list:
    chro_end=chro_sizes[chro]
    # windows is a list of the end of the windows (for this chromosome)
    for window_start in range(0,chro_end,window_size):
        fl=chro+'\t'+str(window_start)+'\t'+str(window_start+window_size)+'\n'
        fakefile.write(fl)
