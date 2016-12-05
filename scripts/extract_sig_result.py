#!/usr/bin/env python
import sys
import os
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests


def main():
	method = sys.argv[1]
	wkdir  = sys.argv[2]
	col_id = int(sys.argv[3])
	col_p  = int(sys.argv[4])
	out = sys.argv[5]
	files = os.listdir(wkdir)
	res = pd.DataFrame()
	for file in files:
		print("Processing", file)
		tb = pd.read_table(os.path.join(wkdir, file), sep='\t', header=0)
		# ignore NA in p-value
		if tb.iloc[:,col_p].isnull().sum() > 0:
			print("WARNING: NA found in P-values")
			tb = tb[~tb.iloc[:,col_p].isnull()]
		# get qvalue
		tb['q'] = multipletests(tb.iloc[:,col_p], method='fdr_bh')[1]
		# q < 0.1
		tb = tb[tb.q<0.1]
		tb['tumor'] = file.split('.')[0]
		tb['method'] = method
		# id p q method tumor
		tb = tb.loc[:, (tb.columns.values[col_id],
			tb.columns.values[col_p], 'q', 'method', 'tumor')]
		tb.columns = ('id', 'p', 'q', 'method', 'tumor')
		print("{} sig. bins".format(tb.shape[0]))
		# append res
		res = res.append(tb, ignore_index=True)
	# return res
	res.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
	main()

