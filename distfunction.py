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
	''' Double Schechter function in logarithmic form.
		Args:
			m: log base 10 of the parameter to bin by (e.g, mass or luminosity)
			mstar1/mstar2: log base 10 of the Schechter parameter normalization
			phistar1/phistar2: log base 10 of the Schechter normalization
			alpha1/alpha2: faint-end Schechter slope
		Output:
			phi: density per logarithmic interval Mpc^{-3} dex^{-1}
	'''
	from numpy import log, exp
	A = pow(10, phistar1) * pow(10, (m - mstar1)*(1 + alpha1)) * exp(- pow(10, m - mstar1))
	B = pow(10, phistar2) * pow(10, (m - mstar2)*(1 + alpha2)) * exp(- pow(10, m - mstar2))
	return log(A + B)

def densityBin(m, vol, bins=51, **kwargs):
	''' Distribution function built by binning data, then returning the bin density and dividing by binwidth to removing binning dependence.
		Args:
			m: parameter to bin by (e.g., mass or luminosity)
			volume: volume of space containing objects of interest (i.e, comoving volume)
		Output:
			phi: density in Mpc^{-3} dex^{-1}
	'''
	from numpy import histogram, diff
	counts, binedges = histogram(m, bins=bins, **kwargs)
	binwidth = diff(binedges)
	return log(counts/vol/binwidth)

