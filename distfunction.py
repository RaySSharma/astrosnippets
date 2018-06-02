#!/usr/bin/env python

def schecterFunc(m, mstar, phistar, alpha):
	''' Schechter (1976) function in logarithmic form. Returns density as a function of m in units of Mpc^{-3} log(m)^{-1}
		Args:
			m: log base 10 of the parameter plotted over (e.g, mass or luminosity)
			mstar: log base 10 of the Schechter parameter normalization
			phistar: log base 10 of the Schechter normalization
			alpha: faint-end Schechter slope 
	'''
	from numpy import log, exp
	return log(10) * pow(10, phistar) * pow(10, (m - mstar)*(1 + alpha)) * exp(- pow(10, m - mstar))

def doubleSchechterFunc(m, mstar1, mstar2, phistar1, phistar2, alpha):
	''' Double Schechter function in logarithmic form. Returns density as a function of m in units of Mpc^{-3} log(m)^{-1}
		Args:
			m: log base 10 of the parameter plotted over (e.g, mass or luminosity)
			mstar1/mstar2: log base 10 of the Schechter parameter normalization
			phistar1/phistar2: log base 10 of the Schechter normalization
			alpha1/alpha2: faint-end Schechter slope 
	'''
	from numpy import log, exp
	A = log(10) * pow(10, phistar1) * pow(10, (m - mstar1)*(1 + alpha1)) * exp(- pow(10, m - mstar1))
	B = log(10) * pow(10, phistar2) * pow(10, (m - mstar2)*(1 + alpha2)) * exp(- pow(10, m - mstar2))
	return A + B