#!/usr/bin/python
# This takes output from refgene and produces a list of EXONIC REGIONS in a bed format, so we can do nice bedmap operations.

import gzip
import sys

#f = gzip.open('/Users/stephanie/ll/data/genes/refgene.txt.gz','r')
f = open('/Users/stephanie/ll/data/genes/Glist_mm9_Stephanie.txt','r')
exonic = open('/Users/stephanie/ll/data/genes/exonic.bed','w')
intronic = open('/Users/stephanie/ll/data/genes/intronic.bed','w')
# headerline
f.readline()

for line in f:
    gene_ID = line.split()[0]
    chro = line.split()[2]
    strand = line.split()[3]
    # this looks a bit confusing, but the information is only used for recording the start of the gene
    if strand=='+':
        # start of the first exon
        tx_start = line.split()[5].strip(',').split(',')[0]
    else:
        # end of the last exon
        tx_start = line.split()[6].strip(',').split(',')[-1]
    tx_end = line.split()[6].strip(',').split(',')[-1]
    exonCount = int(line.split()[4])
    starts = line.split()[5].strip(',')
    ends = line.split()[6].strip(',')
    if not int(ends.split(',')[-1]) == int(tx_end):
        print ends
        print ends[-1]
        print tx_end
        sys.exit('wtf?')
    if len(starts)!=len(ends):
        print gene_ID
        print starts
        print ends
        sys.exit('wtf?')
    exon_pairs = zip(starts.split(','),ends.split(','))
    intron_ends = starts.split(',')[1:]+[tx_end]
    intron_pairs = zip(ends.split(','),intron_ends)
    for i in range(exonCount-1):
        if strand == '+':
            n_exon = i+1
        else:
            n_exon = exonCount - i -1
        exon = exon_pairs[i]
        exonic_newline = chro+'\t'+exon[0]+'\t'+exon[1]+'\t'+gene_ID+'\t'+str(n_exon)+'\t'+strand+'\t'+tx_start+'\n'
        exonic.write(exonic_newline)

        intron = intron_pairs[i]
        intronic_newline = chro+'\t'+intron[0]+'\t'+intron[1]+'\t'+gene_ID+'\t'+str(n_exon)+'\t'+strand+'\t'+tx_start+'\n'
        intronic.write(intronic_newline)
    exon = exon_pairs[exonCount-1]
    if strand=='+':
        n_exon = exonCount
    else:
        n_exon = 1
    exonic_newline = chro+'\t'+exon[0]+'\t'+exon[1]+'\t'+gene_ID+'\t'+str(n_exon)+'\t'+strand+'\t'+tx_start+'\n'
    exonic.write(exonic_newline)
