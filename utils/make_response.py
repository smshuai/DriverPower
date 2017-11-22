#!/usr/bin/env python
"""
Make required response table from mutations, elements and optionally whitelist regions.
All data should be in the BED format.

Usage:
 $ python ./make_response.py mut.bed element.bed [whitelist.bed]

Data requirements:
1. mut.bed
    * 6 columns ['chrom', 'start', 'end', 'ref', 'alt', 'donorID']
    * no header
    * e.g.,
        chr1	99048	99049	T	C	DO52685
        chr1	238650	238651	A	-	DO228371
        chr1	245377	245378	-	TTTC	DO228322

2. element.bed
    * 4 columns ['chrom', 'start', 'end', 'binID']
    * no header
    * e.g.,
        chr1	565900	566050	element_1
        chr1	566760	566910	element_2

3. whitelist.bed [optional]
    * 3 columns ['chrom', 'start', 'end']
    * no header
    * e.g.,
        chr1	10002	10468
        chr1	12753	12782

"""
import sys
import pandas as pd
from pybedtools import BedTool


if __name__ == '__main__':
    if len(sys.argv) == 3:
        muts = BedTool(sys.argv[1])
        bins = BedTool(sys.argv[2])
        all_binIDs = bins.to_dataframe(names=['chrom2', 'start2', 'end2', 'binID'], usecols=['binID'])\
        .binID.unique()
    elif len(sys.argv) == 4:
        # use the whitelisted regions
        muts = BedTool(sys.argv[1])
        bins = BedTool(sys.argv[2])
        all_binIDs = bins.to_dataframe(names=['chrom2', 'start2', 'end2', 'binID'], usecols=['binID'])\
        .binID.unique()
        callable = BedTool(sys.argv[3])
        bins = bins.intersect(callable)  # only keep callable regions
    else:
        sys.stderr.write('Usage: python make_response.py MUTATION ELEMENT [WHITELIST]')
        sys.exit(1)
    mut_bin = muts.intersect(bins, wa=True, wb=True)
    mut_bin = mut_bin.to_dataframe(names=['chrom1', 'start1', 'end1', 'ref', 'alt', 'donorID',
                                          'chrom2', 'start2', 'end2', 'binID'])
    nmut = mut_bin.binID.value_counts()
    nsample = mut_bin.loc[:, ('binID', 'donorID')].drop_duplicates().binID.value_counts()
    N = muts.to_dataframe(names=['chrom1', 'start1', 'end1', 'ref', 'alt', 'donorID'], usecols=['donorID'])\
        .donorID.unique().shape[0]
    bins_df = bins.to_dataframe(names=['chrom2', 'start2', 'end2', 'binID'])
    bins_df['length'] = bins_df['end2'] - bins_df['start2']
    bin_length = bins_df.pivot_table(index='binID', values='length', aggfunc=sum)['length']
    response = pd.DataFrame({'length': bin_length, 'nMut': nmut, 'nSample': nsample})
    response = response.loc[all_binIDs]
    response['N'] = N
    response = response.fillna(0).astype(int)
    response.index.name = 'binID'
    response.to_csv('./response.tsv', sep='\t')