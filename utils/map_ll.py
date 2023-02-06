#!/usr/bin/env python
u"""
map_ll.py
Original IDL program mapll.pro written by Eric Rignot
Adapted by Tyler Sutterley (04/2015)

CALLING SEQUENCE:
	xy = map_ll(alon, alat, HEM='S')
	x = xy['x']
	y = xy['y']

INPUTS:
	alon: longitude
	alat: latitude
OUTPUTS:
	x: horizontal coordinate polar stereographic projection
	y: vertical coordinate polar stereographic projection
OPTIONS:
	HEM: hemisphere ('S' == south, 'N' == north)

NOTES:
	Snyder, 1982, Map Projections used by the U.S. Geological Survey
	Forward formulas for the ellipsoid

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
	Updated 04/2015: added option for Bamber 5km DEM
	Updated 06/2014: updated comments
	Updated 02/2014: minor update to if statements
	Updated 05/2013: converted to python
"""

def map_ll(alon, alat, HEM='S'):
	import numpy as np

	rad_e = 6378137.0#-- WGS84
	ecc2 = 0.00669437999015#-- WGS84
	ecc = np.sqrt(ecc2)

	#-- Standard parallel - latitude with no distortion = +70/-71.
	if HEM in ('N','n'):#-- Northern Hemisphere
		slat = 70.0
		sn = 1.0
		xlam = 45.0
	elif HEM in ('S','s'):#-- Southern Hemisphere
		slat = 71.0
		sn = -1.0
		xlam = 0.0
	elif HEM in ('GDEM'):#-- Bamber Greenland 5km DEM
		slat = 71.0
		sn = 1.0
		xlam = 39.0

	alat = sn*alat
	alon = sn*alon
	rlat = alat

	t1 = np.tan(np.pi/4.0 - rlat*np.pi/(2.0*180.0)) / \
			((1.0 - ecc*np.sin(rlat*np.pi/180.0))/ \
			(1.0 + ecc*np.sin(rlat*np.pi/180.0)))**(ecc/2.0)

	t2 = np.tan(np.pi/4.0 - slat*np.pi/(2.0*180.0)) / \
			((1.0 - ecc*np.sin(slat*np.pi/180.0)) / \
			(1.0 + ecc*np.sin(slat*np.pi/180.0)))**(ecc/2.0)
	#-- m at standard latitude
	cm = np.cos(slat*np.pi/180.0) / \
		np.sqrt(1.0-ecc2*(np.sin(slat*np.pi/180.0)**2.0))
	#-- radius of latitude circle
	rho = rad_e*cm*t1/t2
	#-- polar stereographic x and y
	x = ( rho*sn*np.sin((alon+xlam)*np.pi/180.0))
	y = (-rho*sn*np.cos((alon+xlam)*np.pi/180.0))

	return {'x':x, 'y':y}
