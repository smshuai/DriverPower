#!/usr/bin/env python
import sys

maf  = sys.argv[1]
meta = sys.argv[2]
out  = sys.argv[3]

outf = open(out, 'w')

with open(meta, 'r') as f:
    sid = []
    for line in f:
        sid.append(line.strip().split("\t")[0])

sys.stderr.write('Read {} sample IDs\n'.format(len(sid)))


N = 0
with open(maf, 'r') as f:
    for line in f:
        name = line.strip().split("\t")[12] # tumor barcode
        if name in sid:
            outf.write(line)
            N += 1

outf.close()
sys.stderr.write("Extract {} mutations\n".format(N))


