#!/usr/bin/env python
"""
calc_sig.py
Adapted for python from IDL calc_sig.pro coded by Sean Swenson
Adapted by Tyler Sutterley (05/2013)

Returns the spatial field for a series of spherical harmonics

CALLING SEQUENCE:
	sig_out = calc_sig(clm1,slm1,lon,lat,LMIN=0,LMAX=60)

INPUTS:
	clm: cosine spherical harmonic coefficients
	slm: sine spherical harmonic coefficients
	lon: longitude
	lat: latitude

OPTIONS:
	LMIN: minimum degree of spherical harmonics
	LMAX: maximum spherical harmonic degrees
	PLM: plm coefficients (if computed outside the function)

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

PROGRAM DEPENDENCIES:
	1. plm_mohlenkamp.py: Computes fully-normalized associated Legendre polynomials 
"""

def calc_sig(clm1,slm1,lon,lat,LMIN=0,LMAX=0,PLM=0):
	import numpy as np
	from plm_mohlenkamp import plm_mohlenkamp

	#-- if LMAX is not specified, will use the size of the input harmonics
	if (LMAX == 0):
		LMAX = np.shape(clm1)[0]-1

	phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]#-- Longitude in radians
	th = (90 -np.squeeze(lat))*np.pi/180.0#-- Colatitude in radians

	thmax = len(th)
	phimax= len(phi)
	sig_out=np.zeros((phimax,thmax))

	#---  Calculate fourier coefficients from legendre coefficients  -----! 
	d_cos=np.zeros((LMAX+1,thmax))
	d_sin=np.zeros((LMAX+1,thmax))
	#-- if plms already computed
	if (np.ndim(PLM) == 0):
		PLM = plm_mohlenkamp(LMAX,np.cos(th))
	
	#--- Truncating harmonics to degree and order LMAX
	#--- removing coefficients below LMIN
	clm = np.zeros((LMAX+1,LMAX+1))
	slm = np.zeros((LMAX+1,LMAX+1))
	clm[LMIN:LMAX+1,0:LMAX+1] = clm1[LMIN:LMAX+1,0:LMAX+1]
	slm[LMIN:LMAX+1,0:LMAX+1] = slm1[LMIN:LMAX+1,0:LMAX+1]
	for k in range(0,thmax):
		d_cos[:,k]=np.sum(PLM[:,:,k]*clm[:,:],axis=0)
		d_sin[:,k]=np.sum(PLM[:,:,k]*slm[:,:],axis=0)

	#---  Final signal recovery from fourier coefficients  ---------------!
	m = np.arange(0,LMAX+1)[:,np.newaxis]
	#-- Calculating cos(m*phi) and sin(m*phi)
	ccos = np.cos(np.dot(m,phi))
	ssin = np.sin(np.dot(m,phi))
	sig_out = np.dot(np.transpose(ccos),d_cos) + np.dot(np.transpose(ssin),d_sin)

	return sig_out
