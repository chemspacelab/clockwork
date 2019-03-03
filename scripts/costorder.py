#!/usr/bin/env python
import scipy.special as scs
import itertools as it
import numpy as np

N = 40 # Number of torsions

def configuration_count(resolution):
	if resolution == 0:
		return 1
	return 2**(resolution - 1)

def combinations_in_group(T, R):
	combinations = []
	for case in it.product(range(R+1), repeat=T):
		if R in case:
			combinations.append(case)
	return combinations

def group_cost(combinations):
	cost = 0
	for case in combinations:
		cost += np.prod(map(configuration_count, case))
	return cost

RS = 5
print ('     ' + ' '.join(['  R:%6d' % _ for _ in range(RS)]))
for T in range(1, 5):
	res = []
	for R in range(RS):
		res.append(scs.binom(N, T)*group_cost(combinations_in_group(T, R)))
	print ('T:%3d' % T + ' '.join(['%10d' % _ for _ in res]))
