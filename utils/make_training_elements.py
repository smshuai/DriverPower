"""
Make training elements from test elements.

Usage:
    $ python ./make_training_elements.py element.bed bed12 50 exclude.bed


Data requirements:
1. element.bed
    * bed4 or bed12
    * no header
    * e.g. bed4:
        chr1	565900	566050	element_1
        chr1	566760	566910	element_2

Output:
    training_elements.bed

"""

import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool


if __name__ == '__main__':
    assert len(sys.argv) == 5, 'Usage: python ./make_training_elements.py TEST_ELEMENTS FORMAT FOLD EXCLUDE'
    test = BedTool(sys.argv[1])
    bedFlag = sys.argv[2]
    fold = int(sys.argv[3])
    exclude = BedTool(sys.argv[4])
    if bedFlag == 'bed4':
        train = pd.DataFrame({'chrom': [], 'start': [], 'end': []})
        for i in range(fold):
            tmp = test.shuffle(genome='hg19', seed=i, excl=sys.argv[1]).\
                to_dataframe(names=['chrom', 'start', 'end', 'name'],
                             usecols=['chrom', 'start', 'end'])
            train = pd.concat([train, tmp])
        train = train.astype({'chrom': str, 'start': np.int_, 'end': np.int_}).sort_values(['chrom','start'])
        train = BedTool.from_dataframe(train.loc[:, ('chrom', 'start', 'end')])
        train = train.subtract(exclude, A=True)
        train = train.to_dataframe(names=['chrom', 'start', 'end'])
        train['name'] = train.apply(lambda x: x.chrom + ':' + str(x.start) + '-' + str(x.end), axis=1)
    elif bedFlag == 'bed12':
        first = True
        for i in range(fold):
            tmp = test.shuffle(genome='hg19', seed=i).intersect(test, v=True, split=True).to_dataframe()
            if first:
                train = tmp.copy()
                first = False
            else:
                train = pd.concat([train, tmp])
        # give element new names
        train['name'] = train.apply(lambda x: x.chrom + ':' + str(x.start) + '-' + str(x.end), axis=1)
    else:
        sys.stderr.write('Unknown input BED format: {}. Must be "bed12" or "bed4"\n'.format(bedFlag))
        sys.exit(1)
    train = train.loc[train.chrom != 'chrY']  # remove chrom Y
    train.drop_duplicates(subset='name', inplace=True)  # remove potential duplicates
    train.to_csv('training_elements.bed', sep='\t', index=False, header=False)  # save bed12 or bed4

