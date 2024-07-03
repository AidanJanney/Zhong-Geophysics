#!/usr/bin/env python
# this code will generate:
# - years.reverse: epochs of the dates in ky
# - ice6g_lm.lmax100.txt: ice6g ice height in meter, up to degree and order 100, with the ice stopped at 100 years ago. The normalization is for spherical harmonic as in Jackson, rather than the geodesy normalization. This is for the input of the old fortran code.
# - ocn1.0_xyt.txt: ocean function in the 1x1 spatial domain
# - ocn1.0_lm.txt: ocean function in the spherical harmonic domain
# - ice6g+ocn.lmax100.txt: ice6g ice load + ocean load, all in meter (equivalent ice height)

from os import path
from ncdf_write import ncdf_write
from gen_clms import gen_clms
from calc_sig import calc_sig
import numpy as np
import netCDF4
from os import path
import re
import math
from scipy import interpolate
import pdb

M_PI=np.pi
lmax=100
clm=np.zeros((lmax+1,lmax+1))
slm=np.zeros((lmax+1,lmax+1))

ddir='/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/ice6g_data'
INDEX_FILE='new.index.txt'

of1='years.reverse'
of2='ice6g_lm.lmax100.txt'
of3='ocn6g_xyt.txt'
of4='ocn6g_lm.txt'
of5='ice6g+ocn.lmax100.txt'


op2 = 1
op3 = 1
op4 = 1
if op2:
	f2=open(of2,'w')
if op3:
	f3=open(of3,'w')
if op4:
	f4=open(of4,'w')
if op2 and op4:
	f5=open(of5,'w')

lon0=np.arange(0,361)
lat0=np.arange(-90,91)

input_files=np.loadtxt(path.join(ddir,INDEX_FILE),dtype=np.str)

n_files=len(input_files)
year=np.empty(49)*np.nan

#ice_x = [1,9,27,48]
for i in range(n_files):
#for i in ice_x:
	data=netCDF4.Dataset(path.join(ddir,input_files[i]),'r')
#        pdb.set_trace()
#        print(data)
        lon = np.array(data['lon'])
        lat = np.array(data['lat'])
        ih = np.array(data['stgit']).T
	### isolate Antarctica:
	#antidx = (lat>=-60.)
	#ih[:,antidx] = 0. # remove Antarctica or the rest of the world
	### 

	ice = np.array(data['sftgif']).T
        land = np.array(data['sftlf']).T # tranpose so that (row,col) = (lon,lat)
        data.close()

        #change for the 122 kyr data file(i=0)
        #pdb.set_trace()
        if i==0:
             ih = np.zeros((len(ih),len(ih[0])))
        #     ice = np.zeros((len(ice),len(ice[0])))

	ocn = 0.*land
	idx = (land<50) & (ice<50)
	ocn[idx] = 1.


	# output 1: 'ice_x.txt' and 'ocn_x.txt'
	#pdb.set_trace()
	if 1:
	        size = len(ih)*len(ih[0])
                ih_0 = ih[:,::-1] 
                ih_1 = ih_0.T 
	        ih_2 = np.reshape(ih_1,(size,1))
                fn_ice = 'ice_x'+str(i)+'.txt' 
	        np.savetxt(path.join(ddir,fn_ice),ih_2,fmt='%14.5e')

		#size = len(ih)*len(ih[0])
                ocn_0 = ocn[:,::-1]
                ocn_1 = ocn_0.T
                ocn_2 = np.reshape(ocn_1,(size,1))
                fn_ocn = 'ocn_x'+str(i)+'.txt'
                np.savetxt(path.join(ddir,fn_ocn),ocn_2,fmt='%14.5e')
#        #-- find numerical values within the filename
        date=re.findall("\d+",input_files[i])
	
	if len(date)==4:
		year[i]=float(date[3])
	if len(date)==5:
		year[i]=float(date[3])+0.1*float(date[4])
# output 3 and 4 : save ocean function in spatial and spherical harmonic domains
	if op3:
		f=interpolate.RectBivariateSpline(lon,lat,ocn)
		ocn0 = f(lon0,lat0)
		ocn0[ocn0<0.5]=0.
		ocn0[ocn0>=0.5]=1.
		lat00=lat0[::-1]
		ocn00=ocn0[:,::-1]
		f3.write('t=%d\n'%i)
		for p in  range(len(lat0)):
			for q in range(len(lon0)):
#				f3.write('%f %f %f\n'%(90.-lat00[p],lon0[q],ocn00[q,p]))
				f3.write('%f %f %f\n'%(90-lat0[180-p],lon0[q],ocn0[q,180-p]))
		
	if op4:
		ylms = gen_clms(ocn,lon,lat,lmax,lmax)
		clm = ylms['clm']
		slm = ylms['slm']
		f4.write('t=%d\n'%i)
		for l in range(lmax+1):
			for m in range(l+1):
				if m==0:
					coef=np.sqrt(4*M_PI)
				else:
					coef=np.sqrt(2*M_PI)*(-1.0)**m
				f4.write('%3d %3d %14.5e %14.5e\n'%(l,m,clm[l,m]*coef,-slm[l,m]*coef))
		A00=clm[0,0]*np.sqrt(4*M_PI)
		

	if i>40:
# check if ice height changes in the last couple thousand years.
		print i, np.sum(ih-ih0)
	ih0 = ih.copy()
	if op2:
		ylms = gen_clms(ih,lon,lat,lmax,lmax)
		clm2 = ylms['clm']
		slm2 = ylms['slm']
		f2.write('\n\nt= %d \n'%i)
		for l in range(lmax+1):
			for m in range(l+1):
                                if m==0:
                                        coef=np.sqrt(4*M_PI)
                                else:
                                        coef=np.sqrt(2*M_PI)*(-1.0)**m
                                f2.write(' l=%3d m=%3d %14.5e  %14.5e\n'%(l,m,clm2[l,m]*coef,-slm2[l,m]*coef))
		I00=clm2[0,0]*np.sqrt(4*M_PI)
# add static ocean load to the ice load
#        pdb.set_trace()
	if op4 and op2:
		CON = -I00/A00
                f5.write('\n\nt= %d \n'%i)
                for l in range(lmax+1):
                        for m in range(l+1):
                                if m==0:
                                        coef=np.sqrt(4*M_PI)
                                else:
                                        coef=np.sqrt(2*M_PI)*(-1.0)**m
                                f5.write(' l=%3d m=%3d %14.5e  %14.5e\n'%(l,m,(clm2[l,m]+CON*clm[l,m])*coef,-(slm2[l,m]+CON*slm[l,m])*coef))
		
   
	
if op2:
	f2.close()
if op3:
	f3.close()
if op4:
	f4.close()
if op4 and op2:
	f5.close()

# sort out the order of the files
if 0:
	idx = np.argsort(year)
	newlist=input_files[idx][::-1]
	np.savetxt(path.join(ddir,'new.index.txt'),newlist,fmt='%s')

# output 1: 'years.reverse'
if 0:
        year.sort()
        year = year[::-1]
	np.savetxt(path.join(ddir,of1),year,fmt='%4.1f')

