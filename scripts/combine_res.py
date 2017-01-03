#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test


def main():
    # input
    res_cadd = pd.read_table(sys.argv[1], sep='\t', header=0, index_col='binID')
    res_eigen = pd.read_table(sys.argv[2], sep='\t', header=0, index_col='binID')
    cutoff_cadd = int(sys.argv[3])
    cutoff_eigen = int(sys.argv[4])
    totcg = pd.read_table(sys.argv[5], sep='\t', header=0, index_col='binID')
    N = int(sys.argv[6])
    res_cadd = tune_cutoff(res_cadd, cutoff_cadd, totcg, N)
    res_eigen = tune_cutoff(res_eigen, cutoff_eigen, totcg, N)
    res_cadd.columns = ['Pval_CADD', 'Qval_CADD']
    res_eigen.columns = ['Pval_EIGEN', 'Qval_EIGEN']
    res_comb = res_cadd.join(res_eigen)

def tune_cutoff(res, cutoff, totcg, N):
    threshold = np.percentile(res.fscore, cutoff)
    res['Pval'] = [binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(res.nMutSample, res.Length*N, res.MuAdj)]
    res = totcg.join(res)
    res.Pval.fillna(1, inplace=True)
    res['Qval'] = multipletests(res.Pval, method='fdr_bh')[1]
    res.sort_values('Pval', inplace=True)
    return res.loc[:, ('Pval', 'Qval')]

if __name__ == '__main__':
    main()
