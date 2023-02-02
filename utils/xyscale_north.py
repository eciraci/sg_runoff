#!/usr/bin/env python
u"""
xyscale_north.py
Original IDL program xyscale_south.pro written by Eric Rignot
Adapted by Tyler Sutterley (01/2015)

This function calculates the scaling factor for a polar stereographic
	projection (ie. SSM/I grid) to correct area calculations. The scaling
	factor is defined (from Snyder, 1982, Map Projections used by the U.S.
	Geological Survey) as:

	k = (mc/m)*(t/tc), where:

	m = cos(lat)/sqrt(1 - e2*sin(lat)^2)
	t = tan(Pi/4 - lat/2)/((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2)
	e2 = 0.006693883 is the earth eccentricity (Hughes ellipsoid)
	e2 = 0.00669437999015d0   WGS84
	e = sqrt(e2)
	mc = m at the reference latitude (+70 degrees)
	tc = t at the reference latitude (+70 degrees)

INPUTS:
	lat: Latitude given in radians
		for 70N = +70.0*np.pi/180.0

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY
	Updated 01/2015: forked for Northern Hemisphere
	Updated 06/2014: updated comments
	Written 06/2013
"""

def xyscale_north(lat):
	import numpy as np

	#-- WGS84 ellipsoid flattening
	flat = 1/298.257223563
	#-- square of the eccentricity of the ellipsoid
	#-- ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
	ecc2 = 2.0*flat - flat**2
	#-- eccentricity of the ellipsoid
	ecc = np.sqrt(ecc2)

	#-- calculate ratio at reference latitude
	m70=np.cos(70.0*np.pi/180.0) / \
		np.sqrt(1.0 - ecc2*np.sin(70.0*np.pi/180.0)**2)
	t70=np.tan(np.pi/4.0 - 70.0*np.pi/(2.0*180.0)) / \
		((1.0 - ecc*np.sin(70.0*np.pi/180.0)) / \
		(1.0 + ecc*np.sin(70.0*np.pi/180.0)))**(ecc/2.0)
	m70_t70 = m70/t70

	m = np.cos(lat)/np.sqrt(1.0 - ecc2*np.sin(lat)**2)
	t = np.tan(np.pi/4.0 - lat/2.0)/((1.0-ecc*np.sin(lat)) / \
		(1.0 + ecc*np.sin(lat)))**(ecc/2.0)
	k = m70_t70*t/m
	#-- modification 04/01/02
	scale = (1.0/k/k)

	return scale
