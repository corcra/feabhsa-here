#!/usr/bin/python

import sys

def overlap(line1,line2):
    start1=line1.split()[1]
    start2=line2.split()[1]
    end1=line1.split()[2]
    end2=line2.split()[2]
    chro1=line1.split()[0]
    chro2=line2.split()[0]
    strand1=line1.split()[5]
    strand2=line2.split()[5]
    if not start1 == start2:
#        print 'starts are different'
        return False
    elif not end1 == end2:
#        print 'ends are different'
        return False
    elif not chro1 == chro2:
#        print 'chros are different'
        return False
    elif not strand1 == strand2:
#        print 'strands are different'
        return False
    else:
#        print 'overlap!!!!'
        return True

if len(sys.argv)<2:
    sys.exit("Requires gene list!")

filename=sys.argv[1]
f=open(filename,"r")
prev_line=f.readline()

outfile=open(filename.replace("_with_dupes",""),"w")
for line in f:
    if not overlap(prev_line,line):
        outfile.write(line)
    prev_line=line
