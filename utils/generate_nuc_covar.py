#!/usr/bin/env python
''' Create nucleotide contents features.
'''
import pandas as pd
import numpy as np
import sys



if __name__ == '__main__':
    assert len(sys.argv) == 3, 'Incorrect number of args, usage: generate_nuc_covar.py PATH_CG PATH_OUT'
    path_cg  = sys.argv[1] # coverage table
    path_out = sys.argv[2] # path of output
    cg = pd.read_table(path_cg, header=0, sep='\t', index_col='binID')
    cg.sort_index(inplace=True)
    cg_sum = cg.sum(axis=1)
    # pct of 3mer cg, 32 features, 32 in total
    cg_percent = cg.apply(lambda x: x/cg_sum * 100, axis=0)
    # log(cg+1), 1 feature, 33 in total
    logcg = np.log10(cg_sum+1)
    logcg.name = 'logCG'
    # pct of 2mer cg, 14 features, 47 in total
    ## mutate at 5'
    two_mer1 = ['TA', 'TC', 'TG', 'TT', 'CA', 'CC', 'CG', 'CT']
    two_mer_dict1 = dict()
    for key in two_mer1:
        two_mer_dict1[key] = [i + key for i in ['A', 'C', 'G', 'T']]
    two_mer_cg1 = dict()
    for key in two_mer1:
        two_mer_cg1[key] = cg.loc[:, two_mer_dict1[key]].sum(axis=1)
    two_mer_cg1 = pd.DataFrame(two_mer_cg1)
    two_mer_cg1.columns = two_mer_cg1.columns.values + '5p' # add suffix to columns
    assert np.array_equal(two_mer_cg1.sum(axis=1), cg_sum)
    ## mutate at 3'
    two_mer2 = ['AT', 'CT', 'GT', 'TT', 'AC', 'CC', 'GC', 'TC']
    two_mer_dict2 = dict()
    for key in two_mer2:
        two_mer_dict2[key] = [key + i for i in ['A', 'C', 'G', 'T']]
    two_mer_cg2 = dict()
    for key in two_mer2:
        two_mer_cg2[key] = cg.loc[:, two_mer_dict2[key]].sum(axis=1)
    two_mer_cg2 = pd.DataFrame(two_mer_cg2)
    two_mer_cg2.columns = two_mer_cg2.columns.values + '3p' # add suffix to columns
    assert np.array_equal(two_mer_cg2.sum(axis=1), cg_sum)
    ## remove duplicate columns
    del two_mer_cg2['CC3p'], two_mer_cg2['TT3p']
    two_mer = pd.concat([two_mer_cg1, two_mer_cg2], axis=1)
    ## convert to pct
    two_mer_percent = two_mer.apply(lambda x: x/cg_sum * 100, axis=0)
    two_mer_percent.fillna(0, inplace=True)

    # combine 99 features
    nuc = pd.concat([logcg, cg_percent, two_mer_percent], axis=1)
    nuc.to_csv(path_out, sep='\t')

    sys.exit(0)
