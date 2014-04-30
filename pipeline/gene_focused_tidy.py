# !/usr/bin/python
# The output of the bedmap operation produces quite a mess. We want to know how many gene_start hits and how many gene body hits each gene has. Not so hard, eh?

import sys

bedmap_mess=open('/Users/stephanie/ll/data/genes/analyse_FP/gene_list_with_overlaps.bed','r')
bedmap_mess=open(sys.argv[1],'r')
outfile=open(sys.argv[2],'w')
#errorfile=open('/Users/stephanie/ll/data/genes/analyse_FP/inconsistent.txt','w')

i=0
for line in bedmap_mess:
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
            where_starts.append(o.split()[1])
            where_ends.append(o.split()[2])
    newline = genedata+'\t'+str(totalcounts)+'\t'+str(sum(nGene_start))+'\t'+str(sum(nBody))+'\t'+';'.join(which)+'\t'+';'.join(where_starts)+'\t'+';'.join(where_ends)
    outfile.write(newline+'\n')
#    if not sum(nGene_start) + sum(nBody) +sum(nOutside)== totalcounts:
  #      print 'huh?', nGene_start, nBody, totalcounts, line
#        print i
#        i=i+1
        #errorfile.write(newline+'\t'+'-'.join(map(str,nOutside))+'\n')
#        sys.exit()
