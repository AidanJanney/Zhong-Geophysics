#!/usr/bin/env python
"""
gen_clms.py
Adapted for python from IDL gen_clms.pro coded by Sean Swenson
Adapted by Tyler Sutterley (05/2013)

Computes the spherical harmonic coefficients of a spatial field

CALLING SEQUENCE:
	Ylms = gen_clms(sig_in, lon, lat, lmax, mmax)
	blm = Ylms['blm']
	clm = Ylms['clm']

INPUTS:
	sig_in: input spatial field
	lon: input longitude
	lat: input latitude
	lmax: maximum spherical harmonic degree
	mmax: maximum spherical harmonic order
OUTPUTS:
	clm: cosine spherical harmonic coefficients
	slm: sine spherical harmonic coefficients

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

PROGRAM DEPENDENCIES:
	1. gen_plms.py: Computes the Fourier coefficients of the 
		associated Legendre functions
"""

def gen_clms(sig_in,lon,lat,lmax,mmax):
	from sys import exit as message
	import numpy as np
	from gen_plms import gen_plms
        

	phi = (np.squeeze(lon)*np.pi/180.0)#-- Longitude in radians
	theta = ((90.0 -np.squeeze(lat))*np.pi/180.0)#-- Colatitude in radians

	#--  Determine size of arrays  -----------------------------------
	thmax=len(theta)
	phimax=len(phi)
	nmax=mmax
	ll = np.arange(0,lmax+1)[:, np.newaxis]#-- lmax+1 to include lmax
	mm = np.arange(0,mmax+1)[:, np.newaxis]#-- mmax+1 to include mmax
	nn = np.arange(0,nmax+1)[:, np.newaxis]#-- nmax+1 to include nmax

	#--  Calculate cos and sin coefficients of signal  ----------------
	mp=np.dot(mm,phi[np.newaxis,:])
        #print sig_in.shape, mp.shape
	dft_cos=np.dot(np.cos(mp),sig_in)
	dft_sin=np.dot(np.sin(mp),sig_in)
	#--  Normalise fourier coefficients  ------------------------------
	dft_cos[0,:]=dft_cos[0,:]/phimax
	dft_cos[1:mmax+1,:]=2.0*dft_cos[1:mmax+1,:]/phimax
	dft_sin[0,:]=dft_sin[0,:]/phimax
	dft_sin[1:mmax+1,:]=2.0*dft_sin[1:mmax+1,:]/phimax

	#--  Calculate cos and sin coefficients of theta component  -------
	#--  Because the function is defined on (0,pi)  -------------------
	#--  it can be expanded in just cosine terms.  --------------------
	#--  this routine assumes that 0 and pi are not included 
	theta_cc=np.zeros((mmax+1,nmax+1))
	theta_sc=np.zeros((mmax+1,nmax+1))
	m_even = np.arange(0,mmax+1,2)
	m_odd = np.arange(1,mmax,2)

	if ((theta[0] == 0.0) and (np.abs(np.pi - theta[thmax-1]) < 1e-5)):
		##(theta(thmax) == pi):
		#--  non-endpoints  ---------------------------------------------
		nt=np.dot(nn,theta[1:thmax-1][np.newaxis,:])
		theta_cc[m_even,:]=2.0*np.dot(dft_cos[m_even,1:thmax-1],np.transpose(np.cos(nt)))
		theta_sc[m_even,:]=2.0*np.dot(dft_sin[m_even,1:thmax-1],np.transpose(np.cos(nt)))
		theta_cc[m_odd,:]=2.0*np.dot(dft_cos[m_odd,1:thmax-1],np.transpose(np.sin(nt)))
		theta_sc[m_odd,:]=2.0*np.dot(dft_sin[m_odd,1:thmax-1],np.transpose(np.sin(nt)))

		#--  endpoints  -------------------------------------------------
		theta_cc[m_even,:]=theta_cc[m_even,:] + \
			np.dot((dft_cos[m_even,0]*np.cos(theta[0]) + \
			dft_cos[m_even,thmax-1]*np.cos(theta[thmax-1]))[:,np.newaxis], \
			np.transpose(nn))
		theta_sc[m_even,:]=theta_sc[m_even,:] + \
			np.dot((dft_sin[m_even,0]*np.cos(theta[0]) + \
			dft_sin[m_even,thmax-1]*np.cos(theta[thmax-1]))[:,np.newaxis], \
			np.transpose(nn))
		theta_cc[m_odd,:]=theta_cc[m_odd,:] + \
			np.dot((dft_cos[m_odd,0]*np.sin(theta[0]) + \
			dft_cos[m_odd,thmax-1]*np.sin(theta[thmax-1]))[:,np.newaxis], \
			np.transpose(nn))			
		theta_sc[m_odd,:]=theta_sc[m_odd,:] + \
			np.dot((dft_sin[m_odd,0]*np.sin(theta[0]) + \
			dft_sin[m_odd,thmax-1]*np.sin(theta[thmax-1]))[:,np.newaxis], \
			np.transpose(nn))			

	elif ((theta[0] != 0.0) and (theta[thmax-1] != np.pi)):
		nt=np.dot(nn,theta[np.newaxis,:])
		theta_cc[m_even,:]=2.0*np.dot(dft_cos[m_even,:],np.transpose(np.cos(nt)))
		theta_sc[m_even,:]=2.0*np.dot(dft_sin[m_even,:],np.transpose(np.cos(nt)))
		theta_cc[m_odd,:]=2.0*np.dot(dft_cos[m_odd,:],np.transpose(np.sin(nt)))
		theta_sc[m_odd,:]=2.0*np.dot(dft_sin[m_odd,:],np.transpose(np.sin(nt)))
	else:
		message('Theta coordinates incompatible')


	#--  Normalize theta fourier coefficients  -----------------------------
	theta_cc[:,0]=theta_cc[:,0]/(2.0*thmax)
	theta_cc[:,1:nmax+1]=theta_cc[:,1:nmax+1]/thmax
	theta_sc[:,0]=theta_sc[:,0]/(2.0*thmax)
	theta_sc[:,1:nmax+1]=theta_sc[:,1:nmax+1]/thmax

	#--  Correct normalization for the incomplete coverage of the sphere  -
	delphi=np.abs(phi[1]-phi[0])
	deltheta=np.abs(theta[1]-theta[0])
	norm=phimax*delphi/(2.0*np.pi)*thmax*deltheta/np.pi
	theta_cc=theta_cc*norm
	theta_sc=theta_sc*norm

	#--  Calculate cos and sin coefficients of Legendre functions  ----
	#--  Expand m = even terms in a cosine series  --------------------
	#--  Expand m = odd terms in a sine series  ------------------------
	#--  Both are stride 2  ----------------------------------------------
	plm=gen_plms(lmax,mmax) 

	#--  Sum theta fourier coefficients  ----------------------------    
	#--  temp is the integral of cos(n theta) cos(k theta) dcos(theta)
	#--  over the interval 0 to pi.  n and k must have like parities
	#---  m = even terms  -------------------------------------------
	kn_even=np.zeros((len(m_even),len(m_even)))
	for n in range(0,nmax+2,2):
		kn_even[:,n/2] = 0.5*(1.0/(1.0-m_even-n) +1.0/(1.0+m_even-n) + \
			1.0/(1.0-m_even+n) +1.0/(1.0+m_even+n))

	kn_odd=np.zeros((len(m_odd),len(m_odd)))
	for n in range(1,nmax+1,2):
		kn_odd[:,(n-1)/2] = 0.5*(1.0/(1-m_odd-n) + 1.0/(1+m_odd-n) + \
			1.0/(1-m_odd+n) + 1.0/(1+m_odd+n))

	l_even = np.arange(0,lmax+1,2)
	l_odd = np.arange(1,lmax,2)
	clm = np.zeros((lmax+1,lmax+1))
	slm = np.zeros((lmax+1,lmax+1))
	for m in range(0,nmax+2,2):
		temp=np.dot(np.transpose(np.squeeze(plm[l_even,m,m_even[:,np.newaxis]])),kn_even)
		clm[l_even,m]=np.dot(np.squeeze(theta_cc[m,m_even[:,np.newaxis]]),np.transpose(temp))
		slm[l_even,m]=np.dot(np.squeeze(theta_sc[m,m_even[:,np.newaxis]]),np.transpose(temp))
		temp=np.dot(np.transpose(np.squeeze(plm[l_odd,m,m_odd[:,np.newaxis]])),kn_odd)
		clm[l_odd,m]=np.dot(np.squeeze(theta_cc[m,m_odd[:,np.newaxis]]),np.transpose(temp))
		slm[l_odd,m]=np.dot(np.squeeze(theta_sc[m,m_odd[:,np.newaxis]]),np.transpose(temp))

	#---  m = odd terms  -------------------------------------------
	#---  kn_ matrices defined differently than for m == even  ------
	kn_even=np.zeros((len(m_even),len(m_even)))
	for n in range(0,nmax+2,2):
		kn_even[:,n/2]=0.5*(-1.0/(1-m_even-n) + 1.0/(1.0+m_even-n) + \
			1.0/(1.0-m_even+n) - 1.0/(1.0+m_even+n)) 

	kn_odd=np.zeros((len(m_odd),len(m_odd)))
	for n in range(1,nmax+1,2):
		kn_odd[:,(n-1)/2]=0.5*(-1.0/(1-m_odd-n) + 1.0/(1.0+m_odd-n) + \
			1.0/(1.0-m_odd+n) -1.0/(1.0+m_odd+n))

	l_even = np.arange(2,lmax+1,2)#-- do not in include l=0
	l_odd = np.arange(1,lmax,2)	
	for m in range(1,nmax+1,2):
		temp=np.dot(np.transpose(np.squeeze(plm[l_even,m,m_even[:,np.newaxis]])),kn_even)
		clm[l_even,m]=np.dot(np.squeeze(theta_cc[m,m_even[:,np.newaxis]]),np.transpose(temp))
		slm[l_even,m]=np.dot(np.squeeze(theta_sc[m,m_even[:,np.newaxis]]),np.transpose(temp))
		temp=np.dot(np.transpose(np.squeeze(plm[l_odd,m,m_odd[:,np.newaxis]])),kn_odd)
		clm[l_odd,m]=np.dot(np.squeeze(theta_cc[m,m_odd[:,np.newaxis]]),np.transpose(temp))
		slm[l_odd,m]=np.dot(np.squeeze(theta_sc[m,m_odd[:,np.newaxis]]),np.transpose(temp))


	#---  Divide by Plm normalization  ------------------------------
	clm[:,0]=clm[:,0]/2.0
	clm[:,1:lmax+1]=clm[:,1:lmax+1]/4.0
	slm[:,0]=slm[:,0]/2.0
	slm[:,1:lmax+1]=slm[:,1:lmax+1]/4.0	
	
	return {'clm':clm, 'slm':slm}
