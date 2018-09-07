""" Prepare data in required format for DriverPower
Input: mutations, elements, whitelisted regions
Output: response table
"""

import logging, sys
import pandas as pd
import numpy as np
from pybedtools import BedTool

logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s | %(levelname)s: %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP-prepare')


def main(path_to_mut, path_to_ele, path_to_callable, path_to_out):
    mut = BedTool(path_to_mut)  # chr1    1230447 1230448 G       A       DO46416
    ele = BedTool(path_to_ele).sort()  # chr10   100002554       100002997       bin1
    whitelist = BedTool(path_to_callable).sort()  # chr1    10002   10468
    N = len(pd.read_table(path_to_mut, header=None, usecols=(5,))[5].unique())  # Get total number of samples
    # Only keep whitelisted elements
    old_cov = ele.total_coverage()
    ele = ele.intersect(whitelist).sort()
    new_cov = ele.total_coverage()
    logger.info("Whitelisted elements: {}/{} bp ({:.2f}%)".format(new_cov, old_cov, new_cov/old_cov*100))
    # Elements intersect mutations
    cnames = ['chrom', 'start', 'end', 'ref', 'alt', 'donor', 'chrom2', 'start2', 'end2', 'binID']
    mut_ele = mut.intersect(ele, wa=True, wb=True).to_dataframe(names=cnames)
    # nMut and nSample per table
    response_tab = mut_ele.pivot_table(index='binID', values='donor', aggfunc=[len, lambda x: len(x.unique())])
    response_tab.columns = ['nMut', 'nSample']
    # Get effective length
    ele_df = ele.to_dataframe(names=['chrom', 'start', 'end', 'binID'])
    ele_df['length'] = ele_df['end'] - ele_df['start']
    eff_length = ele_df.pivot_table(index='binID', values='length', aggfunc=sum)
    response_tab = eff_length.join(response_tab)
    response_tab = response_tab.fillna(0).astype(np.int_)
    response_tab['N'] = N
    response_tab.to_csv(path_to_out, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
