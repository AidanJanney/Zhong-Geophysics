#!/usr/bin/env python
# This code reads the Stokes coefficients of radial (UP) and horizontal (East and North) displacement in meters
# Generate the RSL history for a given location

from os import path
from ncdf_write import ncdf_write
from plm_mohlenkamp import plm_mohlenkamp
from gen_clms import gen_clms
from calc_sig import calc_sig
from calc_hori import calc_hori
import numpy as np
from os import path
import re

lmax=100
M_PI = np.pi
fact=1.e-2 # convert cm to meters

nepochs = 122
lon = np.array([360-59.53, 360-59.53])
lat = np.array([13.07, 13.07])

clm = np.zeros((lmax+1,lmax+1))
slm = np.zeros((lmax+1,lmax+1))
uclm = np.zeros((lmax+1,lmax+1))
uslm = np.zeros((lmax+1,lmax+1))
plm = plm_mohlenkamp(lmax,np.sin(lat*np.pi/180.0))

rsl = np.zeros(nepochs)*np.nan

for count in range(nepochs):
	infile='PW.ice6g.'+str(101+count)
	data = np.loadtxt(infile,dtype=np.float)
	# data format is (l, m, clm, slm, uclm, uslm, hclm, hslm), i.e. degree, order, geoid stokes coeff, uplift Stokes, and horizontal Stokes

	i=0
	for l in range(lmax+1):
		for m in range(l+1):
			if (l!=data[i,0]) or (m!=data[i,1]):
				print "something is wrong!", l, m, data[i,0], data[i,1]
				stop
			if m==0:
				coef = 1./np.sqrt(4*M_PI)
			else:
				coef = 1./(np.sqrt(2*M_PI)*(-1)**m)
			clm[l,m] = data[i,2]*coef
			slm[l,m] = -data[i,3]*coef
			uclm[l,m] = data[i,4]*coef
			uslm[l,m] = -data[i,5]*coef
			i = i+1

	grate = calc_sig(clm,slm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
	urate = calc_sig(uclm,uslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
	rsl[count] = grate[0,0] - urate[0,0]
print grate.shape

print "Done calculation..."

# save text files:
if 1:
	otxt1='RSL_N_U_incom_Barbo.txt'
	np.savetxt(otxt1,rsl,fmt='%14.5e  ')

