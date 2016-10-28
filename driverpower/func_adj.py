''' Functional adjustment module for DriverPower
'''

import tabix
import logging


# create  logger
logger = logging.getLogger('FUNC ADJ')
logger.setLevel(logging.INFO)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s',
    datefmt='%m/%d/%Y %H:%M:%S')
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(ch)


def get_eigen(coding=True):
	''' Get eigen scores with tabix
	'''
	pass


def get_cadd():
	''' Get CADD scores with tabix
	'''
	pass


def func_adj(method='eigen', cutoff=0.8):
	''' Main wrapper for functional adjustment
	'''
	support_method = ['eigen']
	assert method in support_method, 'Invalid functional score method. Must be chosen from {}'.format(support_method)
	if method == 'eigen':
		pass