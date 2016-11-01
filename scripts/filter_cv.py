#!/usr/bin/env python
from driverpower import load
from driverpower import preprocess
import sys

if __name__ == '__main__':
    path_ct = sys.argv[1]
    path_coverage = sys.argv[2]
    path_cv = sys.argv[3]
    path_out = sys.argv[4]
    ct = load.load_count(path_ct)
    cg = load.load_coverage(path_coverage)
    keep, recur = preprocess.get_filter(ct, cg, return_recur=True)
    out = open(path_out, 'w')
    with open(path_cv, 'r') as f:
        out.write(f.readline()) # header
        for line in f:
            arr = line.strip().split("\t")
            binID = arr[0]
            if recur.loc[binID] >= 2:
                out.write(line)
    out.close()
