#!/usr/bin/env python
''' Expand a bed, like from genes to exons
Usage:
    python expand_bed.py < in.bed > out.bed
'''
import fileinput
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) # avoid broken pipe error

def process_line(line):
    gene = line.strip().split('\t')
    chrom = gene[0]
    binID = gene[3]
    start = gene[11].split(',')
    start.pop() # remove last element
    start = [int(s) + int(gene[1]) for s in start]
    length = gene[10].split(',')
    length.pop()
    end = [int(s) + int(l) for s, l in zip(start, length)]
    uniq_binID = [ binID + '::' + str(i) for i in list(range(len(start)))]
    exons = zip([chrom]*len(start), start, end, [binID]*len(start), uniq_binID) # chrom start end binID uniq_binID
    assert len(start) == int(gene[9])
    return exons


for line in fileinput.input():
    exons = process_line(line)
    exons = ['\t'.join(str(i) for i in e) + '\n' for e in exons]
    sys.stdout.write(''.join(exons))
