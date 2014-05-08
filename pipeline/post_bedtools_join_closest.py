# !/usr/bin/python
# Bedtools spits out ties on new lines. This script puts them back on the same line. >:|

import sys

firstline = sys.stdin.readline()
current_ID = firstline.split()[0]
distance_buffer=firstline.split()[2]
gene_buffer=[firstline.split()[1]]
strand_buffer={firstline.split()[3]}
for line in sys.stdin:
    ID=line.split()[0]
    distance=line.split()[2]
    gene=line.split()[1]
    strand=line.split()[3]
    if ID==current_ID:
        gene_buffer.append(gene)
        strand_buffer.add(strand)
    else:
        buffered_line=';'.join(gene_buffer)+'\t'+distance_buffer+'\t'+';'.join(strand_buffer)
        print buffered_line
        gene_buffer=[gene]
        current_ID = ID
        distance_buffer = distance
        strand_buffer = {strand}

buffered_line=';'.join(gene_buffer)+'\t'+distance_buffer+'\t'+strand
print buffered_line
