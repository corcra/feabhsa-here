# !/usr/bin/python
# The output of the bedmap operation produces quite a mess. We want to know how many gene_start hits and how many gene body hits each gene has. Not so hard, eh?

import sys

i=0
for line in sys.stdin:
    nGene_start=[]
    nBody=[]
    nOutside=[]
    which = []
    where_starts = []
    where_ends = []
    genedata=line.split('|')[0]
    totalcounts=int(line.split('|')[1])
    overlaps = line.split('|')[2:]
    for o in overlaps:
        if len(o)>1:
            nGene_start.append(int(o.split()[12]))
            nBody.append(1*(o.split()[13]!='0'))
            nOutside.append(int(o.split()[14]))
            which.append(o.split()[3])
            # only record the location if it's in the gene body...
            if o.split()[13]!='0':
                where_starts.append(o.split()[1])
                where_ends.append(o.split()[2])
    outline = genedata+'\t'+str(totalcounts)+'\t'+str(sum(nGene_start))+'\t'+str(sum(nBody))+'\t'+';'.join(which)+'\t'+';'.join(where_starts)+'\t'+';'.join(where_ends)
    print outline
