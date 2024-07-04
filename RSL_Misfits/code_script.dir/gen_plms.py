#!/usr/bin/env python
"""
gen_plms.py
Adapted for python from IDL gen_plms.pro coded by Sean Swenson
Adapted by Tyler Sutterley (04/2013)

Computes the Fourier coefficients of the associated Legendre
	functions needed by gen_clms.py

CALLING SEQUENCE:
	plm = gen_plms(lmax,mmax)

INPUTS:
	lmax: maximum spherical harmonic degree
	mmax: maximum spherical harmonic order
OUTPUTS:
	plm: Fourier coefficients

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)
"""

def gen_plms(lmax,mmax):
	import numpy as np

	plm = np.zeros((lmax+1,lmax+1,lmax+1))
	l_even = np.arange(0,lmax+1,2)
	l_odd = np.arange(1,lmax,2)
	m_even = np.arange(0,mmax+1,2)
	m_odd = np.arange(1,mmax,2)

	#--  First compute m=0, m=1 terms  -------------------------------
	#---  Compute m = 0, l = even terms  ------------------------------
	plm[l_even,0,0] = 1.0
	plm[l_even,0,2] = (l_even*(l_even+1))*plm[l_even,0,0]/(l_even*(l_even+1)-2)
	for j in range(2,lmax,2):#-- equivalent to 2:lmax-2
		plm[l_even,0,j+2] = (2.0*(l_even*(l_even+1)-j**2)*plm[l_even,0,j] + \
			((j-2)*(j-1)-l_even*(l_even+1))*plm[l_even,0,j-2]) / \
			(l_even*(l_even+1)-(j+2)*(j+1))


	#---  Special case for j = 0 fourier coefficient  --------------------
	plm[l_even,0,0] = plm[l_even,0,0]/2.0
	#---  Normalize overall sum to 2 for m eq 0  ----------------------
	norm = np.zeros((len(l_even)))
	for j in range(0,lmax+2,2):#-- equivalent to 0:lmax
		norm[l_even/2] = norm[l_even/2]+plm[l_even,0,j] * \
			np.dot(np.squeeze(plm[l_even[:, np.newaxis],0,m_even]), \
			(1.0/(1-j-m_even)+1.0/(1+j-m_even) + \
			1.0/(1-j+m_even) +1.0/(1+j+m_even)))/2.0

	norm = np.sqrt(norm/2.0)
	for l in range(0,lmax+2,2):#-- equivalent to 0:lmax
		plm[l,0,:] = plm[l,0,:]/norm[l/2]


	#---  Compute m = 0, l = odd terms  -------------------------------
	plm[l_odd,0,1] = 1.0
	plm[l_odd,0,3] = (2-l_odd*(l_odd+1))*plm[l_odd,0,1]/(6-l_odd*(l_odd+1))
	for j in range(3,lmax-1,2):#-- equivalent to 3:lmax-3
		plm[l_odd,0,j+2] = (2.0*(l_odd*(l_odd+1)-j**2)*plm[l_odd,0,j]+((j-2)*(j-1) - \
			l_odd*(l_odd+1))*plm[l_odd,0,j-2])/(l_odd*(l_odd+1)-(j+2)*(j+1))


	#---  Normalize overall sum to 2 for m == 0  ----------------------
	norm = np.zeros((len(l_odd)))
	for j in range(1,lmax+1,2):#-- equivalent to 1:lmax-1
		norm[(l_odd-1)/2] = norm[(l_odd-1)/2]+plm[l_odd,0,j] * \
			np.dot(np.squeeze(plm[l_odd[:, np.newaxis],0,m_odd]), \
			(1.0/(1-j-m_odd)+1.0/(1+j-m_odd) + \
			1.0/(1-j+m_odd) +1.0/(1+j+m_odd)))/2.0

	norm = np.sqrt(norm/2.0)
	for l in range(1,lmax+1,2):#-- equivalent to 1:lmax-1
		plm[l,0,:] = plm[l,0,:]/norm[(l-1)/2]
	

	#---  Compute m = 1, l = even terms  ------------------------------
	plm[l_even,1,0] = 0.0
	plm[l_even,1,2] = 1.0
	for j in range(2,lmax,2):#-- equivalent to 2:lmax-2
		plm[l_even,1,j+2] = (2.0*(l_even*(l_even+1)-j**2-2)*plm[l_even,1,j] + \
			((j-2)*(j-1)-l_even*(l_even+1))*plm[l_even,1,j-2]) / \
			(l_even*(l_even+1)-(j+2)*(j+1))

	#---  Normalize overall sum to 4 for m == 1  ----------------------
	#---  this is different than that for the cosine series  --------------
	norm = np.zeros((len(l_even)))
	for j in range(0,lmax+2,2):#-- equivalent to 0:lmax
		norm[l_even/2] = norm[l_even/2]+plm[l_even,1,j] * \
			np.dot(np.squeeze(plm[l_even[:, np.newaxis],1,m_even]), \
			(-1.0/(1-j-m_even)+1.0/(1+j-m_even) + \
			1.0/(1-j+m_even) - 1.0/(1+j+m_even)))/2.0

	norm = np.sqrt(norm/4.0)
	for l in range(0,lmax+2,2):#-- equivalent to 0:lmax
		plm[l,1,:] = plm[l,1,:]/norm[l/2]

	#---  Compute m = 1, l = odd terms  -------------------------------
	plm[l_odd,1,1] = 1.0
	plm[l_odd,1,3] = 3.0*(l_odd*(l_odd+1)-2)*plm[l_odd,1,1]/(l_odd*(l_odd+1)-6)
	for j in range(3,lmax-1,2):#-- equivalent to 3:lmax-3
		plm[l_odd,1,j+2] = (2.0*(l_odd*(l_odd+1)-j**2-2)*plm[l_odd,1,j] + \
			((j-2)*(j-1)-l_odd*(l_odd+1))*plm[l_odd,1,j-2]) / \
			(l_odd*(l_odd+1)-(j+2)*(j+1))

	#---  Normalize overall sum to 4 for m == 1  ----------------------
	norm = np.zeros((len(l_odd)))
	for j in range(1,lmax+1,2):#-- equivalent to 1:lmax-1
		norm[(l_odd-1)/2] = norm[(l_odd-1)/2] + plm[l_odd,1,j] * \
			np.dot(np.squeeze(plm[l_odd[:, np.newaxis],1,m_odd]), \
			(-1.0/(1-j-m_odd)+1.0/(1+j-m_odd) + \
			1.0/(1-j+m_odd) - 1.0/(1+j+m_odd)))/2.0

	norm = np.sqrt(norm/4.0)
	for l in range(1,lmax+1,2):#-- equivalent to 1:lmax-1
		plm[l,1,:] = plm[l,1,:]/norm[(l-1)/2]
	

	#---  Compute coefficients for m > 0  -------------------------------
	#---  Long integers must be used for loop indeces  ----------------
	#---  m = 0 terms on rhs have different normalization  -------------
	m = 0
	#---  m = 0, l = even terms  -----------------------------------------
	for l in range(m,lmax-1):#-- equivalent to m:lmax-2
		plm[l+2,m+2,m_even] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1)/ \
			np.float(2*l+1))*plm[l,m,m_even] - \
			np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
			np.float(2*l+5))*plm[l+2,m,m_even] + \
			np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)/2.0) * \
			plm[l,m+2,m_even]) / np.sqrt(np.float(l+m+4) * \
			np.float(l+m+3)/np.float(2*l+5)/2.0)
	
	#---  m = 0, l = odd terms  ---------------------------------------
	for l in range(m+1,lmax-1):#-- equivalent to m+1:lmax-2
		plm[l+2,m+2,m_odd] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1)/ \
			np.float(2*l+1))*plm[l,m,m_odd] - \
			np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
			np.float(2*l+5))*plm[l+2,m,m_odd] + \
			np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)/2.0) * \
			plm[l,m+2,m_odd]) / np.sqrt(np.float(l+m+4) * \
			np.float(l+m+3)/np.float(2*l+5)/2.0)

	#---  m = even, > 2, l = even terms  ------------------------------
	for m in range(2,lmax,2):#-- equivalent to 2:lmax-2 
		for l in range(m,lmax,2):#-- equivalent to m:lmax-2 
			plm[l+2,m+2,m_even] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1) / \
				np.float(2*l+1))*plm[l,m,m_even] - \
				np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
				np.float(2*l+5))*plm[l+2,m,m_even] + \
				np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)) * \
				plm[l,m+2,m_even]) / np.sqrt(np.float(l+m+4) * \
				np.float(l+m+3)/np.float(2*l+5))

		#---  m = even, > 2, l = odd terms  -------------------------------
		for l in range(m+1,lmax-1,2):
			plm[l+2,m+2,m_odd] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1) / \
				np.float(2*l+1))*plm[l,m,m_odd] - \
				np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
				np.float(2*l+5))*plm[l+2,m,m_odd] + \
				np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)) * \
				plm[l,m+2,m_odd]) / np.sqrt(np.float(l+m+4) * \
				np.float(l+m+3)/np.float(2*l+5))


	#---  m = odd, > 1, l = even terms  -------------------------------
	for m in range(1,lmax-1,2):#-- equivalent to 1:lmax-3
		for l in range(m+1,lmax-1,2):#-- equivalent to m+1,lmax-2
			plm[l+2,m+2,m_even] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1) / \
				np.float(2*l+1))*plm[l,m,m_even] - \
				np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
				np.float(2*l+5))*plm[l+2,m,m_even] + \
				np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)) * \
				plm[l,m+2,m_even]) / np.sqrt(np.float(l+m+4) * \
				np.float(l+m+3)/np.float(2*l+5))

	#---  m = odd, > 1, l = odd terms  --------------------------------
		for l in range(m,lmax-1,2):#-- equivalent to m:lmax-2
			plm[l+2,m+2,m_odd] = (np.sqrt(np.float(l+m+2)*np.float(l+m+1) /
				np.float(2*l+1))*plm[l,m,m_odd] - \
				np.sqrt(np.float(l-m+1)*np.float(l-m+2) / \
				np.float(2*l+5))*plm[l+2,m,m_odd] + \
				np.sqrt(np.float(l-m)*np.float(l-m-1)/np.float(2*l+1)) * \
				plm[l,m+2,m_odd]) / np.sqrt(np.float(l+m+4) * \
				np.float(l+m+3)/np.float(2*l+5))

	return plm
