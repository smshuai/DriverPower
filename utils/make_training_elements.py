"""
Make training elements from test elements.

Usage:
    $ python ./make_training_elements.py element.bed 50 exclude.bed


Data requirements:
1. element.bed
    * 4 columns ['chrom', 'start', 'end', 'binID']
    * no header
    * e.g.,
        chr1	565900	566050	element_1
        chr1	566760	566910	element_2

Output:
    training_elements.tsv

"""

import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool


if __name__ == '__main__':
    assert len(sys.argv) == 4, 'Usage: python ./make_training_elements.py TEST_ELEMENTS FOLD EXCLUDE'
    test = BedTool(sys.argv[1])
    fold = int(sys.argv[2])
    exclude = BedTool(sys.argv[3])
    train = pd.DataFrame({'chrom': [], 'start': [], 'end': []})
    for i in range(fold):
        tmp = test.shuffle(genome='hg19', seed=i, excl=sys.argv[1]).\
            to_dataframe(names=['chrom', 'start', 'end', 'binID'],
                         usecols=['chrom', 'start', 'end'])
        train = pd.concat([train, tmp])
    train = train.astype({'chrom': str, 'start': np.int_, 'end': np.int_}).sort_values(['chrom','start'])
    train = BedTool.from_dataframe(train.loc[:, ('chrom', 'start', 'end')])
    train = train.subtract(exclude, A=True)
    train = train.to_dataframe(names=['chrom', 'start', 'end'])
    train['binID'] = train.apply(lambda x: x.chrom + ':' + str(x.start) + '-' + str(x.end), axis=1)
    train = train.loc[train.chrom!='chrY']
    train.drop_duplicates(inplace=True)
    train.to_csv('training_elements.bed', sep='\t', index=False, header=False)