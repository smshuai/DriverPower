#!/usr/bin/env python
''' This script is used in beds2cv.sh. Collapse cv table by binID.
Args:
    cvIN  - sys.argv[1]
    cvOUT - sys.argv[2]
Input cv table (tsv with header) format:
    Column 1: binID
    Column 2-ncol: one cv per column
Output cv table format:
    Column 1: binID
    Column 2-ncol: one cv per column
Note:
    nrow(input) >= nrow(output)
'''
import pandas as pd
import sys

if __name__ == '__main__':
    assert len(sys.argv) == 3, 'Incorrect number of args. Usage pivot_cv.py vIN cvOUT'
    # read
    cv = pd.read_csv(sys.argv[1], header=0, sep='\t')
    cv.pivot_table(index='binID', aggfunc=sum).to_csv(sys.argv[2], header=True, sep='\t', index=True)
    sys.exit(0)
    

