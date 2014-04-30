# !/usr/bin/python
# Bedtools spits out ties on new lines. This script puts them back on the same line. >:|

import sys

infile_path = sys.argv[1]
outfile_path = sys.argv[2]

distances_file=open(infile_path,'r')
outfile=open(outfile_path,'w')

firstline = distances_file.readline()
current_ID = firstline.split()[0]
distance_buffer=firstline.split()[2]
gene_buffer=[firstline.split()[1]]
for line in distances_file:
    ID=line.split()[0]
    distance=line.split()[2]
    gene=line.split()[1]
    strand=line.split()[3]
    if ID==current_ID:
        gene_buffer.append(gene)
    else:
        buffered_line=';'.join(gene_buffer)+'\t'+distance_buffer+'\t'+strand+'\n'
        outfile.write(buffered_line)
        gene_buffer=[gene]
        current_ID = ID
        distance_buffer = distance

buffered_line=';'.join(gene_buffer)+'\t'+distance_buffer+'\t'+strand+'\n'
outfile.write(buffered_line)
