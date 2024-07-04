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
from plm_mohlenkamp import plm_mohlenkamp
import numpy as np
import netCDF4
from os import path
import re
import math
from scipy import interpolate
import pdb

M_PI=np.pi
lmax=100

ddir_PW='./density_4400/'

clm = np.zeros((lmax+1,lmax+1))
slm = np.zeros((lmax+1,lmax+1))
uclm = np.zeros((lmax+1,lmax+1))
uslm = np.zeros((lmax+1,lmax+1))
rsl = np.zeros(122)*np.nan

ddir='/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/ice6g_data_48'
INDEX_FILE='new.index.txt'



lon0=np.arange(0,361)
lat0=np.arange(-90,91)

input_files=np.loadtxt(path.join(ddir,INDEX_FILE),dtype=np.str)
n_files=len(input_files)

ddir1='/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/ice6g_data_122'
year_122=np.loadtxt(path.join(ddir1,'years_122.reverse'),dtype=np.str)
year_48=np.loadtxt(path.join(ddir1,'years_48.reverse'),dtype=np.str)

#read longer history data
fn='ICE-6G_C_IceThickness_1deg.nc'
input_file=path.join(ddir1,fn)
data=netCDF4.Dataset(input_file,'r')
print(data)
lon = np.array(data['Lon'])
lat = np.array(data['Lat'])
ih_total = np.array(data['stgit'])
time = np.array(data['Time'])


for count in range(len(time)):
        
	#find the correct index for ocn from 26 ka history
	diff = np.zeros(len(input_files))
	for i_ocn in range(len(input_files)):
        	diff[i_ocn] = abs(float(year_122[count])-float(year_48[i_ocn]))
        #pdb.set_trace()
	idx_ocn = np.where(diff==diff.min())[0]
 	print(count,idx_ocn[0])

	#read data of 26ka
	data=netCDF4.Dataset(path.join(ddir,input_files[idx_ocn[0]]),'r')
        lon = np.array(data['lon'])
        lat = np.array(data['lat'])
        #ih = np.array(data['stgit']).T
        ice = np.array(data['sftgif']).T
        land = np.array(data['sftlf']).T # tranpose so that (row,col) = (lon,lat)
        topo = np.array(data['Topo']).T
        data.close()
        ocn = 0.*land
        idx = (land<50) & (ice<50)
        ocn[idx] = 1.

        # calucate [O]_00
	ylms = gen_clms(ocn,lon,lat,lmax,lmax)
        A00 = ylms['clm'][0,0]

        # calculate [I]_00
        ih = ih_total[count]
	ylms = gen_clms(ih.T,lon,lat,lmax,lmax)
        I00 = ylms['clm'][0,0]

        # calculate [(N-U)*O]_00
        infile='PW.ice6g.'+str(101+count)
        data = np.loadtxt(path.join(ddir_PW,infile),dtype=np.float)
        # data format is (l, m, clm, slm, uclm, uslm, hclm, hslm), i.e. degree, order, geoid stokes coeff, uplift Stokes, and horizontal Stokes
        i = 0
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
                        i = i + 1

        plm = plm_mohlenkamp(lmax,np.sin(lat*np.pi/180.0))
        geoid = calc_sig(clm,slm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*1.e-2
        uplift = calc_sig(uclm,uslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*1.e-2

        ylms = gen_clms((geoid-uplift)*ocn,lon,lat,lmax,lmax)
        NUO00 = ylms['clm'][0,0]

        rsl[count] = (-893./1027.*I00 - NUO00)/A00


# svae text files
if 1:
        otxt1 = 'RSL_c.txt'
        np.savetxt(path.join(ddir,otxt1),rsl,fmt='%14.5e  ')


