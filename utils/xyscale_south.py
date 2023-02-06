#!/usr/bin/env python
u"""
xyscale_south.py
Original IDL program xyscale_south.pro written by Eric Rignot
Adapted by Tyler Sutterley (06/2013)

This function calculates the scaling factor for a polar stereographic
	projection (ie. SSM/I grid) to correct area calculations. The scaling
	factor is defined (from Snyder, 1982, Map Projections used by the U.S.
	Geological Survey) as:

	k = (mc/m)*(t/tc), where:

	m = cos(lat)/sqrt(1 - e2*sin(lat)^2)
	t = tan(Pi/4 - lat/2)/((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2)
	e2 = 0.006693883 is the earth eccentricity (Hughes ellipsoid)
	e2 = 0.00669437999015d0   (WGS84)
	e = sqrt(e2)
	mc = m at the reference latitude (-71 degrees)
	tc = t at the reference latitude (-71 degrees)

	The ratio mc/tc is precalculated and stored in the variable m71_t71.

INPUTS:
	lat: South latitude given in radians (i.e. negative degrees)
		for -71N = +71.0*np.pi/180.0

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY
	Updated 06/2014: updated comments
	Written 06/2013
"""

def xyscale_south(lat):
	import numpy as np

	m71_t71 = 1.9390295644e0
	# m71_t71 = 1.9390300458e0
	#-- square of the eccentricity of the ellipsoid
	#-- ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
	#flat = 1/298.257223563#-- WGS84 ellipsoid flattening
	#ecc2 = 2.0*flat - flat**2
	ecc2= 0.00669437999015e0#-- modification Aug. 26, 2000.
	#-- eccentricity of the ellipsoid
	ecc = np.sqrt(ecc2)

	#m71=np.cos(71.0*np.pi/180.0) / \
	#	np.sqrt(1.0 - ecc2*np.sin(71.0*np.pi/180.0)**2)
	#t71=np.tan(np.pi/4.0 - 71.0*np.pi/(2.0*180.0)) / \
	#	((1.0 - ecc*np.sin(71.0*np.pi/180.0)) / \
	#	(1.0 + ecc*np.sin(71.0*np.pi/180.0)))**(ecc/2.0)
	#m71_t71 = m71/t71

	m = np.cos(lat)/np.sqrt(1.0 - ecc2*np.sin(lat)**2)
	t = np.tan(np.pi/4.0 - lat/2.0)/((1.0-ecc*np.sin(lat)) / \
		(1.0+ecc*np.sin(lat)))**(ecc/2.0)
	k = m71_t71*t/m
	scale = (1.0/k/k)#-- modification 04/01/02

	return scale
