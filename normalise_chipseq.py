#!/usr/bin/python
# Normalising modifications on H3 (ChIPseq data) - K4me1, K4me3, K27ac
# Takes bedGraph files!

import gzip

# Files, options etc.
bedgraph_folder='/Users/stephanie/ll/data/ChIPseq/bedgraphs_strandsep'
chro_data='/Users/stephanie/ll/data/mm9.chromSizes'
window_size=1000

chro_sizes_dat=open(chro_data).readlines()
chro_sizes_list=[(chro.split()[0],int(chro.split()[1])) for chro in chro_sizes_dat]
chro_sizes={chro:size for (chro,size) in chro_sizes_list}

# this is the order chromosomes appear in sorted bed files...
chro_list=['chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr8','chr9','chrM','chrX','chrY']

def normalise(window_counts,window_size,infile,outfile):


def window_overlap_count(start,end,count,window_start,window_end):
    


def get_window_counts(reffile,window_size):

    for chro in chro_list:
        first_line=reffile.readline()
        file_chro=fist_line.split()[0]
        chro_end=chro_sizes[chro]
        # windows is a list of the end of the windows (for this chromosome)
        windows=range(window_size,chro_end,window_size)
        if not windows[-1] == chro_end:
            # the last window is a bit smaller, ok
            windows.append(chro_end)
        for window in windows:
            if filechro==chro:
                while

    # window_counts is a list of tuples, of the form (chro,start(of window),end(of window),counts)
    window_counts=[]

    # start it off
    start=int(first_line.split()[1])
    end=int(first_line.split()[2])
    # define the boundaries of the window
    window_start=start
    window_end=start+window_size
    window_chro=chro
    window_count=0.0

    # progress...
    i=0
    for line in reffile:
        if i%50000==0:
            print i
        i=i+1

        chro=line.split()[0]
        start=int(line.split()[1])
        end=int(line.split()[2])
        count=float(line.split()[3])
        if window_chro==chro:
        # still on the same chromosome
            if not overlap(start,end,window_start):
                # make a new window
                # record previous window info
                window_counts.append((window_chro,window_start,window_end,window_count))
                # here's the new window (yep windows are continuous)
                window_start=window_end+1
                if not window_start+1000 > chro_sizes[window_chro]:
                    window_end = window_start+1000
                else:
                    # at the end of the chromosome (so the window is smaller...)
                    window_end = chro_sizes[window_chro]
            else:
                # overlap between current region and the window
                window_count = window_count + window_overlap_count(start,end,count,window_start,window_end)
        # different chromosome: don't need to check overlap... save the old window, make a new window on a new chromosome
        else:
            window_counts.append((window_chro,window_start,window_end,window_count))
            window_start=0
            window_end=1000
            window_chro=chro
            window_count = window_overlap_count(start,end,count,window_start,window_end)
    # final line/end

total_H3=open('')
