#!/bin/pytohn

import sys

for line in sys.stdin:
    bypipe = line.split('|')
    # first line is kinda the general everything-keeper, it will contain the union of all reigons (e.g. the final region)... the purpose is to enable the bedmap operation to detect overlaps
    overall = bypipe[0]
    overall_chr = overall.split()[0]
    overall_start = int(overall.split()[1])
    overall_end = int(overall.split()[2])
    overall_ID = overall.split()[3]
    for time_overlap in bypipe[1:]:
        start = int(time_overlap.split()[1])
        end = int(time_overlap.split()[2])
        if start < overall_start:
            overall_start = start
        if end > overall_end:
            overall_end = end
   new_overall = overall_chr+'\t'+str(overall_start)+'\t'+str(overall_end)+'\t'+overall_ID
   newline = '|'.join([new_overall]+bypipe[1:])
   print newline


