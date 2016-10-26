#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np


if len(sys.argv) != 5:
    sys.stderr.write('Incorrect number of args. Usage pivot_cg.py PATH_TO_CG BED OUT_CG OUT_TOTCG')
    sys.exit(1)

# read
cg = pd.read_csv(sys.argv[1], header=0, sep='\t', dtype={'chrom':str})
del cg['chrom'], cg['start'], cg['end']
bed = pd.read_table(sys.argv[2], header=None) # 4th column should be binID

# Collapse by binID
res = cg.pivot_table(index='binID', aggfunc=sum)

# missing bins
binID = bed[3].unique()
nmiss = binID.shape[0] - res.shape[0]

if nmiss > 0:
    sys.stderr.write('INFO | {} bins have no covered bases. Filling with zeros\n'.format(nmiss))
    res = res.loc[binID]
    res.fillna(0, inplace=True)

assert np.array_equal(res.index, binID), 'ERROR | binID in cg table does not match binID in BED'

res = res.astype(np.int)
res.sort_index(inplace=True)
res.to_csv(sys.argv[3], header=True, sep='\t')

# Calculate total cg (effective size)
totcg = res.sum(1)
totcg.name = 'totcg'
totcg.to_csv(sys.argv[4], sep='\t', header=True)

sys.stderr.write('INFO | Total number of bases covered = {}\n'.format(res.sum().sum()))
sys.stderr.write('INFO | Mean number of bases covered in bins = {}\n'.format(int(np.mean(res.sum(1)))))
sys.stderr.write('INFO | SD number of bases covered in bins = {}\n'.format(int(np.std(res.sum(1)))))
