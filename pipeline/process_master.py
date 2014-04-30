#!/usr/bin/python

import sys
import gzip
import re

#base=sys.argv[1]+'/'
#datatype=sys.argv[2]
infile_path=sys.argv[1]
outfile_path=sys.argv[2]
datatype=sys.argv[3]

print 'process_master.py: Processing', infile_path, 'into', outfile_path

# what to expect in the master file?
if datatype=='FP':
    times=['0','2','5','12.5','25','50']
elif datatype=='TRP':
    times=['DMSO','12.5','25','50']

n=0
#master=gzip.open(base+datatype+'_dREG_regions.bed.gz','r')
infile=gzip.open(infile_path,'r')
outfile=gzip.open(outfile_path,'w')

header='chr\tstart\tend\tdREG_id\t'+'\t'.join(times)+'\tinflation\twhen_smallest'+'\n'
outfile.write(header)
for line in infile:
    time_pattern=[]
    n=n+1
    chro=line.split()[0]
    uniq_ID='dREG_'+datatype+'_'+str(n)
    uniq_start=line.split()[1]
    uniq_end=line.split()[2]
    # gonna look for the smallest region which contributed to the final region
    smallest = int(uniq_end)-int(uniq_start)
    when_smallest=re.sub('_\d+[|]+[chr]*[1-9MXY]*[1-9]*','',line.split()[3])
#    when_smallest=re.sub('[|]+[chr]+[1-9MXY][1-9]*','',line.split()[3])
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
                when_smallest=re.sub('_\d$','',segment.split()[3])
        else:
            time_pattern.append('0')
    # compared to the smallest hit in this region, how much has our final region expanded?
    final_size=int(uniq_end)-int(uniq_start)
    inflation=float(final_size)/smallest
    uniq_line=chro+'\t'+uniq_start+'\t'+uniq_end+'\t'+uniq_ID+'\t'+'\t'.join(time_pattern)+'\t'+str(inflation)+'\t'+when_smallest+'\n'
    outfile.write(uniq_line)
