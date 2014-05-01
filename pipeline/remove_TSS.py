#!/usr/bin/python

import sys

gb=open(sys.argv[1],'r')
gene_body=open('gene_body.bed','w')

for line in gb:
    bybar=line.split('|')
    chro=bybar[0].split()[0]
    gene_start=bybar[0].split()[1]
    gene_end=bybar[0].split()[2]
    gene_ID=bybar[0].split()[3]
    score=bybar[0].split()[4]
    strand=bybar[0].split()[5]
    line_fragment=chro+'\t'+gene_start+'\t'
    for cut in bybar[1:]:
        if len(cut)>1:
            cut_start=cut.split()[1]
            cut_end=cut.split()[2]
            cut_ID=cut.split()[3]
            if int(cut_start)<int(gene_start):
                # basically just moving the start of the gene up...
                line_fragment=chro+'\t'+cut_end+'\t'
            else:
                line_fragment=line_fragment+cut_start+'\t'+gene_ID+'\t'+score+'\t'+strand+'\t'+cut_ID+'\n'+chro+'\t'+cut_end+'\t'
    final_line=line_fragment+gene_end+'\t'+gene_ID+'\t'+score+'\t'+strand+'\t'+'NA'+'\n'
    gene_body.write(final_line)
