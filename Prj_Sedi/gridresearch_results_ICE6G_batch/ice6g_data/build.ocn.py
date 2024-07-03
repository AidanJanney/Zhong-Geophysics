#!/usr/bin/env python
# this code will generate:
# - ocn1.0_xyt.txt: ocean function in the 1x1 spatial domain
# - ocn1.0_lm.txt: ocean function in the spherical harmonic domain

from os import path
from tyler_library.ncdf_write import ncdf_write
from grace_library.gen_clms import gen_clms
from grace_library.calc_sig import calc_sig
import numpy as np
import netCDF4
from os import path
import re
import math
from scipy import interpolate

M_PI=np.pi
lmax=256

ddir='/Users/geruo/data.dir/ICE6G.dir'
INDEX_FILE='21ka.1ka.index.txt'

lon0=np.arange(0,361)
lat0=np.arange(-90,91)

of3='nocn6g_xyt.txt'
of4='nocn6g_lm.txt'

op3 = 1
op4 = 1
if op3:
        f3=open(of3,'w')
if op4:
        f4=open(of4,'w')

input_files=np.loadtxt(path.join(ddir,INDEX_FILE),dtype=np.str)
n_files=len(input_files)
year=np.empty(48)*np.nan

for i in range(n_files):
	data=netCDF4.Dataset(input_files[i],'r')
        lon = np.array(data['lon'])
        lat = np.array(data['lat'])
        ih = np.array(data['stgit']).T
	ice = np.array(data['sftgif']).T
        land = np.array(data['sftlf']).T # tranpose so that (row,col) = (lon,lat)
        data.close()
	ocn = 0.*land
	idx = (land<50) & (ice<50)
	ocn[idx] = 1.
#        #-- find numerical values within the filename
#        date=re.findall("\d+",input_files[i])
#	if len(date)==5:
#		year[i]=float(date[4])
#	if len(date)==6:
#		year[i]=float(date[4])+0.1*float(date[5])
# output 3 and 4 : save ocean function in spatial and spherical harmonic domains
	if op3:
		f=interpolate.RectBivariateSpline(lon,lat,ocn)
		ocn0 = f(lon0,lat0)
		ocn0[ocn0<0.5]=0.
		ocn0[ocn0>=0.5]=1.
		lat00=lat0[::-1]
		ocn00=ocn0[:,::-1]
#		f3.write('t=%d\n'%i)
		for p in  range(len(lat0)):
			for q in range(len(lon0)):
#				f3.write('%f %f %f\n'%(90.-lat00[p],lon0[q],ocn00[q,p]))
				f3.write('%f %f %f\n'%(90-lat0[180-p],lon0[q],ocn0[q,180-p]))
		
	if op4:
		ylms = gen_clms(ocn,lon,lat,lmax,lmax)
		clm = ylms['clm']
		slm = ylms['slm']
#		f4.write('t=%d\n'%i)
		for l in range(lmax+1):
			for m in range(l+1):
				if m==0:
					coef=np.sqrt(4*M_PI)
				else:
					coef=np.sqrt(2*M_PI)*(-1.0)**m
				f4.write('%3d %3d %14.5e %14.5e\n'%(l,m,clm[l,m]*coef,-slm[l,m]*coef))
		A00=clm[0,0]*np.sqrt(4*M_PI)
	
if op3:
	f3.close()
if op4:
	f4.close()

