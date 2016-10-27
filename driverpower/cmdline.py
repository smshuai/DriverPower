"""
"""
from driverpower.load import load_all
from driverpower.preprocess import preprocess


def main():
	cg_test, ct_test, X_test, cg_train, ct_train, X_train = load_all()
	preprocess()


if __name__ == '__main__':
	main()