#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np


if len(sys.argv) != 3:
    sys.stderr.write('Incorrect number of args. Usage pivot_cg.py PATH_TO_CG OUT')
    sys.exit(1)

# read
cg = pd.read_csv(sys.argv[1], header=0, sep='\t', dtype={'chrom':str})
del cg['chrom'], cg['start'], cg['end']
res = cg.pivot_table(index='binID', aggfunc=sum)
res.to_csv(sys.argv[2], header=True, sep='\t')

sys.stderr.write('INFO: Total number of bases covered = {}\n'.format(res.sum().sum()))
sys.stderr.write('INFO: Mean number of bases covered in bins = {}\n'.format(int(np.mean(res.sum(1)))))
sys.stderr.write('INFO: SD number of bases covered in bins = {}\n'.format(int(np.std(res.sum(1)))))
