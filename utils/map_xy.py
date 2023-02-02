#!/usr/bin/env python
u"""
map_xy.py
Original IDL program mapxy.pro written by Eric Rignot
Adapted by Tyler Sutterley (04/2015)

CALLING SEQUENCE:
	ll = map_xy(x, y, HEM='S')
	lon = ll['lon']
	lat = ll['lat']

INPUTS:
	x: horizontal coordinate polar stereographic projection
	y: vertical coordinate polar stereographic projection
OUTPUTS:
	alon: longitude
	alat: latitude
OPTIONS:
	HEM: hemisphere ('S' == south, 'N' == north)

NOTES:
	Snyder (1982) Map Projections used by the U.S. Geological Survey
		Inverse formulas for the Ellipsoid
	Adams (1921) Latitude developments connected with geodesy and
		cartography with tables

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
	Updated 04/2015: added option for Bamber 5km DEM
	Updated 06/2014: updated comments for references and expanded
		Adams series solution to order 4
	Updated 02/2014: minor update to if statements
	Updated 06/2013: converted to python
"""

def map_xy(x,y, HEM='S'):
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

	#-- size of the input array
	sz = np.ndim(x)
	x=x*1.0
	y=y*1.0

	#-- finding if points are exactly on the pole
	if (sz == 0):
		if ((x == 0.0) and (y == 0.0)):
			#-- making y slightly off pole to prevent division error
			y = 1e-20
	else:
		#-- number of points exactly on the pole
		count = np.count_nonzero((x == 0.0) & (y == 0.0))
		#-- making y slightly off pole to prevent division error
		if (count > 0):
			y[np.nonzero((x == 0) & (y == 0.0))] = 1e-20

	#-- distance rho
	rho=np.sqrt(x**2 + y**2)
	#-- m at standard latitude
	cm = np.cos(slat*np.pi/180.0) / \
		np.sqrt(1.0 - ecc2*(np.sin(slat*np.pi/180.0)**2))
	t = np.tan(np.pi/4.0 - slat*np.pi/(2.0*180.0)) / \
		((1.0 - ecc*np.sin(slat*np.pi/180.0)) / \
		(1.0 + ecc*np.sin(slat*np.pi/180.0)))**(ecc/2.0)
	t = rho*t/(rad_e*cm)
	#-- conformal latitude
	chi = (np.pi/2.0)-2.0*np.arctan(t)
	#-- inverse formula for calculating lat in terms in chi
	#-- Snyder equation 3-5 and Adams (1921) Page 85 phi-chi
	#-- series solution to avoid iterating to convergence
	alat = chi + \
		((ecc2/2.) + (5.*ecc2**2./24.) + (ecc2**3./12.) + (13.*ecc2**4./360.))*np.sin(2.*chi) + \
		((7.*ecc2**2./48.) + (29.*ecc2**3./240.) + (811.*ecc2**4./11520.))*np.sin(4.*chi) + \
		((7.*ecc2**3./120.) + (81.*ecc2**4./1120.))*np.sin(6.*chi) + \
		(4279.*ecc2**4./161280.)*np.sin(8.*chi)
	alat = (sn*alat)*(180.0/np.pi)
	xpr=sn*x
	ypr=sn*y
	#-- arctan2 computes the arctangent of xpr/ypr with range (-pi,pi)
	#-- it will output the correct quadrant and is equivalent to:
	# alon = 2.0*np.arctan(xpr/(np.sqrt(xpr**2+ypr**2)+xpr))*180.0/np.pi - sn*xlam
	alon = np.arctan2(xpr,-ypr)*(180.0/np.pi) - sn*xlam
	alon = sn*alon

	# Fixing longitudes to be 0:360
	#-- input data is points
	if (sz == 0):
		if (alon < 0):
			alon = alon+360.
		elif (alon > 360.0):
			alon = alon-360.
	else: #-- input data is arrays
		ind = np.nonzero(alon < 0.0)
		count = np.count_nonzero(alon < 0.0)
		if (count != 0):
			alon[ind]=alon[ind]+360.0
		ind = np.nonzero(alon > 360.0)
		count = np.count_nonzero(alon > 360.0)
		if (count != 0):
			alon[ind]=alon[ind]-360.0

	return {'lon':alon,'lat':alat}
