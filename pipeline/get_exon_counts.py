#!/bin/python

import sys

for line in sys.stdin:
    geneid = line.split()[3]
    exon_hits = line.strip().split(';')[1:]
    exon_count = 0
    for hit in exon_hits:
        if hit == geneid:
            exon_count = exon_count + 1
    outline = geneid+"\t"+str(exon_count)
    print outline
