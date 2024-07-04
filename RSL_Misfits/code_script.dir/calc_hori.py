#!/usr/bin/env python
"""
calc_hori.py
Adapted from Tyler Sutterley (05/2013) by Geruo A

Returns the horizontal displacements for a series of spherical harmonics

CALLING SEQUENCE:
	east, north = calc_hori(clm1,slm1,lon,lat,LMIN=0,LMAX=60)

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

def calc_hori(clm1,slm1,lon,lat,LMIN=0,LMAX=0,PLM=0):
	import numpy as np
	from plm_mohlenkamp import plm_mohlenkamp

	#-- if LMAX is not specified, will use the size of the input harmonics
	if (LMAX == 0):
		LMAX = np.shape(clm1)[0]-1

	phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]#-- Longitude in radians
	th = (90 -np.squeeze(lat))*np.pi/180.0#-- Colatitude in radians

	thmax = len(th)
	phimax= len(phi)
	east=np.zeros((phimax,thmax))
	north=np.zeros((phimax,thmax))

	#---  Calculate fourier coefficients from legendre coefficients  -----! 
	d_cos1=np.zeros((LMAX+1,thmax))
	d_sin1=np.zeros((LMAX+1,thmax))
        d_cos2=np.zeros((LMAX+1,thmax))
        d_sin2=np.zeros((LMAX+1,thmax))

	#-- if plms already computed
	if (np.ndim(PLM) == 0):
		PLM = plm_mohlenkamp(LMAX,np.cos(th))

	FLM1 = np.zeros((LMAX+1,LMAX+1,thmax))
	FLM2 = np.zeros((LMAX+1,LMAX+1,thmax))
	for m in range(1,LMAX+1):
		for k in range(thmax):
			if (th[k] == 0.) or (th[k]==np.pi):
#				th[k]=1.e-6
				continue
			FLM1[:,m,k] = m*PLM[:,m,k]/np.sin(th[k])
	for l in range(1,LMAX+1):
		for m in range(l+1):
			for k in range(thmax):
				if (th[k] ==0.) or (th[k]==np.pi):
#					th[k]=1.e-6
					continue
				FLM2[l,m,k] = -(l*np.cos(th[k])*PLM[l,m,k] - np.sqrt((2.*l+1.)/(2.*l-1.)*(l*l-m*m))*PLM[l-1,m,k])/np.sin(th[k])

	#--- Truncating harmonics to degree and order LMAX
	#--- removing coefficients below LMIN
	clm = np.zeros((LMAX+1,LMAX+1))
	slm = np.zeros((LMAX+1,LMAX+1))
	clm[LMIN:LMAX+1,0:LMAX+1] = clm1[LMIN:LMAX+1,0:LMAX+1]
	slm[LMIN:LMAX+1,0:LMAX+1] = slm1[LMIN:LMAX+1,0:LMAX+1]
	for k in range(0,thmax):
		d_cos1[:,k]=np.sum(FLM1[:,:,k]*slm[:,:],axis=0)
		d_sin1[:,k]=-np.sum(FLM1[:,:,k]*clm[:,:],axis=0)
                d_cos2[:,k]=np.sum(FLM2[:,:,k]*clm[:,:],axis=0)
                d_sin2[:,k]=np.sum(FLM2[:,:,k]*slm[:,:],axis=0)

	#---  Final signal recovery from fourier coefficients  ---------------!
	m = np.arange(0,LMAX+1)[:,np.newaxis]
	#-- Calculating cos(m*phi) and sin(m*phi)
	ccos = np.cos(np.dot(m,phi))
	ssin = np.sin(np.dot(m,phi))
	east = np.dot(np.transpose(ccos),d_cos1) + np.dot(np.transpose(ssin),d_sin1)
	north = np.dot(np.transpose(ccos),d_cos2) + np.dot(np.transpose(ssin),d_sin2)

	return east, north
