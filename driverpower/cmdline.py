"""
"""
from driverpower.load import load_all
from driverpower.preprocess import preprocess
from driverpower.feature_select import fselect

def main():
	# load all the data. fnames is feature names
	cg_test, ct_test, X_test, cg_train, ct_train, X_train, fnames = load_all()
	
	# get response, filter and scale.
	## gnames is filtered binIDs in test data
	## grecur is filtered recur in test data (index by gnames)
	X_train, ybinom_train, X_test, ybinom_test, gnames, grecur = preprocess(cg_test, ct_test, X_test, cg_train,
    ct_train, X_train, len_threshold, recur_threshold, scaler_type)
	
	# feature selection
	X_train, X_test, fscores = fselect(X_train, X_test, ybinom_train, fnames, method='rndlasso')

	# model


if __name__ == '__main__':
	main()