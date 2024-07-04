#!/usr/bin/env python

from os import path
from ncdf_write import ncdf_write
from plm_mohlenkamp import plm_mohlenkamp
from gen_clms import gen_clms
from calc_sig import calc_sig
import numpy as np
import netCDF4
from os import path
import re
import math
import pdb
from scipy import interpolate

M_PI=np.pi
lmax=100
clm=np.zeros((lmax+1,lmax+1))
slm=np.zeros((lmax+1,lmax+1))
re=6371*1.e3  #6371*1.e3

ddir='/home/khuan/gridresearch_results_ICE6G_batch/ice6g_data'
INDEX_FILE='new.index.txt'

input_files=np.loadtxt(path.join(ddir,INDEX_FILE),dtype=np.str)
n_files=len(input_files)
year=np.empty(n_files)*np.nan

clm = np.zeros((lmax+1,lmax+1))
slm = np.zeros((lmax+1,lmax+1))
uclm = np.zeros((lmax+1,lmax+1))
uslm = np.zeros((lmax+1,lmax+1))
rsl = np.zeros(n_files)*np.nan
rsl_eustatic = np.zeros(n_files)*np.nan
ice_volume = np.zeros(n_files)
for count in range(n_files):
#for count in [49,50]:
	data=netCDF4.Dataset(path.join(ddir,input_files[count]),'r')
        lon = np.array(data['lon'])
        lat = np.array(data['lat'])
        ih = np.array(data['stgit']).T
	ice = np.array(data['sftgif']).T
        land = np.array(data['sftlf']).T # tranpose so that (row,col) = (lon,lat)
	topo = np.array(data['Topo']).T
        data.close()

        if count==0:
             ih = np.zeros((len(ih),len(ih[0])))

	ocn = 0.*land
	idx = (land<50) & (ice<50)
	ocn[idx] = 1.

        if count == 0:
                ocn_LGM = ocn.copy()

        if count == 1:
                #print(ocn_LGM)
                ocn = ocn_LGM    # set LGM ocean function at t=-116


	# calucate [O]_00
	ylms = gen_clms(ocn,lon,lat,lmax,lmax)
	A00 = ylms['clm'][0,0]

	# calculate [I]_00
	ylms = gen_clms(ih,lon,lat,lmax,lmax)
	I00 = ylms['clm'][0,0]

	# calculate [(N-U)*O]_00
	ddir1 = '.'#'../check_code/case_incompr_Lambeck1/'#'../check_code/CASE001/'
        infile='PW.ice6g.'+str(101+count)
        data = np.loadtxt(path.join(ddir1,infile),dtype=np.float)
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
	
	#pdb.set_trace()
	rsl[count] = (-893./1027.*I00 - NUO00)/A00
        ice_volume[count] = re**2*I00*np.sqrt(4*M_PI) 

        rsl_eustatic[count] = (-893./1027.*I00)/A00
# svae text files
if 1:
	otxt1 = 'RSL_c_vm5a_compr.txt'
	np.savetxt(otxt1,rsl,fmt='%14.5e  ')

if 1:
        otxt3 = 'RSL_eustatic_vm5a_compr.txt'
        np.savetxt(otxt3,rsl_eustatic,fmt='%14.5e  ')



if 1:
        otxt2 = 'ice_volume.txt'
        np.savetxt(otxt2,ice_volume,fmt='%14.5e  ')
