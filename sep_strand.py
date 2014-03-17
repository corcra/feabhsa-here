#!/usr/bin/python

import sys

if len(sys.argv)<2:
    sys.exit('give me a bed file, please')

filename = sys.argv[1]

plusfile = open(filename+'.plus','w')
minusfile = open(filename+'.minus','w')

i=0
for line in open(filename,'r'):
    if i%5000==0:
        print i
    i = i+1
    strand = line.split()[5]
    if strand=='+':
        plusfile.write(line)
    else:
        minusfile.write(line)
