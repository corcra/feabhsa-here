#!/usr/bin/python

import sys
import gzip


datatype=sys.argv[1]
base='/Users/stephanie/ll/results/'+datatype+'/'

# this currently has no puprose
gene_list_path='/Users/stephanie/ll/data/genes/gene_list.bed'
gene_list=open(gene_list_path,'r').readlines()
key_vals = [(line.split()[3],(int(line.split()[1]),int(line.split()[2]),line.split()[5])) for line in gene_list]
gene_dict={ID:val for (ID,val) in key_vals}

# what to expect in the master file?
if datatype=='FP':
    times=['0','2','5','12.5','25','50']
elif datatype=='TRP':
    times=['DMSO','12.5','25','50']

n=0
master=gzip.open(base+'dREG_regions.bed.gz','r')
new_master=gzip.open(base+'dREG_regions_uniq.bed.gz','w')

header='chr\tstart\tend\tdREG_id\t'+'\t'.join(times)+'\tinflation'+'\n'
new_master.write(header)
for line in master:
    time_pattern=[]
    n=n+1
    chro=line.split()[0]
    uniq_ID='dREG_'+datatype+'_'+str(n)
    uniq_start=line.split()[1]
    uniq_end=line.split()[2]
    # gonna look for the smallest region which contributed to the final region
    smallest = int(uniq_end)-int(uniq_start)
    for time in times:
        if datatype+'_'+time+'_' in line:
            # yep, this region appeared at this timepoint... mark it down!
            time_pattern.append('1')
            # split the line into timepoint-data (by |), ask which one has the timepoint we want, then isolate that segment
            where=[s.find(datatype+'_'+time+'_')>0 for s in line.split('|')].index(True)
            segment=line.split('|')[where]
            start=segment.split()[1]
            # expand the boundaries
            if int(start)<int(uniq_start):
                uniq_start = start
            end=segment.split()[2]
            if int(end)>int(uniq_end):
                uniq_end = end
            size=int(end)-int(start)
            if size<smallest:
                smallest = size
        else:
            time_pattern.append('0')
    # compared to the smallest hit in this region, how much has our final region expanded?
    final_size=int(uniq_end)-int(uniq_start)
    inflation=float(final_size)/smallest
    uniq_line=chro+'\t'+uniq_start+'\t'+uniq_end+'\t'+uniq_ID+'\t'+'\t'.join(time_pattern)+'\t'+str(inflation)+'\n'
    new_master.write(uniq_line)
