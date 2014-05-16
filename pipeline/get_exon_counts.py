#!/bin/python

import sys

for line in sys.stdin:
    geneid = line.split()[3]
    exon_hits = line.strip().split(';')[1:]
    exon_count = 0
    length = int(line.split()[2])-int(line.split()[1])
    for hit in exon_hits:
        if hit == geneid:
            exon_count = exon_count + 1
    exon_density = float(exon_count)*1000/length
    outline = geneid+"\t"+str(exon_count+"\t"+str(exon_density)
    print outline
