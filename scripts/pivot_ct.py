#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv) != 3:
    sys.stderr.write('Incorrect number of args. Usage pivot_ct.py PATH_TO_CT OUT')

# read
ct = pd.read_csv(sys.argv[1], header=None, sep='\t')
# col12 non-unique binID, col7 donor id, col 8 categ. 
ct = ct.pivot_table(index=[12, 7, 8], values=[0], aggfunc=len).reset_index()
ct.columns = ['binID', 'sid', 'categ', 'ct']
ct.to_csv(sys.argv[2], header=True, index=False, sep='\t')

