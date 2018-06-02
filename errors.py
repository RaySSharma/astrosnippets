#!/usr/bin/env python

def partial_derivative(func, x, var=0, point=[]):
	'''
	Returns the partial derivative of func wrt the nth variable at the given point, but not with respect to x.
	'''
	from scipy.misc import derivative
	args = point[:]
	def wraps(a):
	    args[var] = a
	    return func(x, *args)
	return derivative(wraps, point[var], dx = 1e-8)

def sigma(func, x, params, variances):
	'''
	Returns the total 1-sigma errors of a function given parameter values and their errors. Calculates the square root of the sum of the squared partial derivatives times the respective variances.
	'''
	from numpy import sum, sqrt
	num_var = range(len(params))
	arg = list( pow(partial_derivative(func, x, var=i, point=params),2) * variances[i] for i in num_var )
	return sqrt(sum(arg, axis=0))