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

lon0=np.arange(0,361)
lat0=np.arange(-90,91)

ddir_ice6g='../ice6g_data_48'
input_files_ice6g=np.loadtxt(path.join(ddir_ice6g,'new.index.txt'),dtype=np.str)
n_files_ice6g=len(input_files_ice6g)

ddir='../W12/for_Shijie'
INDEX_FILE='ice_fn1.txt'

input_files0=np.loadtxt(path.join(ddir,INDEX_FILE),dtype=np.str)
input_files = input_files0[::-1]
n_files=len(input_files)
nepoch=len(input_files)

clm = np.zeros((lmax+1,lmax+1))
slm = np.zeros((lmax+1,lmax+1))
uclm = np.zeros((lmax+1,lmax+1))
uslm = np.zeros((lmax+1,lmax+1))
rsl = np.zeros(nepoch)*np.nan

#directory of PW.ice6g.*
ddir_PW='../check_code/case_incompr_W12_b'

year=np.loadtxt(path.join(ddir,'years.reverse'),dtype=np.str)
year_ice6g=np.loadtxt(path.join(ddir_ice6g,'years_48.reverse'),dtype=np.str)

#read longer history data
#fn='ICE-6G_C_IceThickness_1deg.nc'
#input_file=path.join(ddir1,fn)
#data=netCDF4.Dataset(input_file,'r')
#print(data)
#lon = np.array(data['Lon'])
#lat = np.array(data['Lat'])
#ih_total = np.array(data['stgit'])
#time = np.array(data['Time'])


for count in range(len(input_files)):
        
	#find the correct index for ocn from 26 ka history
	diff = np.zeros(len(input_files_ice6g))
	for i_ocn in range(len(input_files_ice6g)):
        	diff[i_ocn] = abs(float(year[count])-float(year_ice6g[i_ocn]))
        #pdb.set_trace()
	idx_ocn = np.where(diff==diff.min())[0]
 	print(count,idx_ocn[0])

	#read data of 26ka
	data_ice6g=netCDF4.Dataset(path.join(ddir_ice6g,input_files_ice6g[idx_ocn[0]]),'r')
        lon_ice6g = np.array(data_ice6g['lon'])
        lat_ice6g = np.array(data_ice6g['lat'])
        #ih = np.array(data['stgit']).T
        ice_ice6g = np.array(data_ice6g['sftgif']).T
        land_ice6g = np.array(data_ice6g['sftlf']).T # tranpose so that (row,col) = (lon,lat)
        topo_ice6g = np.array(data_ice6g['Topo']).T
        data_ice6g.close()
        ocn = 0.*land_ice6g
        idx = (land_ice6g<50) & (ice_ice6g<50)
        ocn[idx] = 1.

        # calucate [O]_00
	ylms = gen_clms(ocn,lon_ice6g,lat_ice6g,lmax,lmax)
        A00 = ylms['clm'][0,0]

        #ih = ih_total[count]
        data_ih = np.loadtxt(path.join(ddir,input_files[count]))
        
	column1 = data_ih[:,0].reshape((256,512))
        lon = column1[0,:]

        column2 = data_ih[:,1].reshape((256,512))
        lat = column2[:,0]

        ih = data_ih[:,2].reshape((len(lat),len(lon))).T
	
        # calculate [I]_00
	ylms = gen_clms(ih,lon,lat,lmax,lmax)
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

        plm = plm_mohlenkamp(lmax,np.sin(lat_ice6g*np.pi/180.0))
        geoid = calc_sig(clm,slm,lon_ice6g,lat_ice6g,LMIN=0,LMAX=lmax,PLM=plm)*1.e-2
        uplift = calc_sig(uclm,uslm,lon_ice6g,lat_ice6g,LMIN=0,LMAX=lmax,PLM=plm)*1.e-2

        ylms = gen_clms((geoid-uplift)*ocn,lon_ice6g,lat_ice6g,lmax,lmax)
        NUO00 = ylms['clm'][0,0]

        rsl[count] = (-893./1027.*I00 - NUO00)/A00


# svae text files
if 1:
        otxt1 = 'RSL_c.txt'
        np.savetxt(path.join(ddir_PW,otxt1),rsl,fmt='%14.5e  ')


