#!/usr/bin/env python
# this code reads the Stokes coefficients of radial (UP) and horizontal (East and North) rates in mm/yr

from os import path
#from tyler_library.ncdf_write import ncdf_write
#from grace_library.plm_mohlenkamp import plm_mohlenkamp
#from grace_library.gen_clms import gen_clms
#from grace_library.calc_sig import calc_sig
#from grace_library.calc_hori import calc_hori
from plm_mohlenkamp import plm_mohlenkamp
from ncdf_write import ncdf_write
from gen_clms import gen_clms
from calc_sig import calc_sig
from calc_hori import calc_hori
import numpy as np
from os import path
import re
import pdb

lmax=32
M_PI = np.pi
fact=1.e-2

casetime = '345'#'49'
caseid0 = '60t'#'11a'
#ty = 'tps'  #uplift
#ty = 'pttl' #geoid
#ddir = '../check_code/CASE001/'
ddir = '../check_code/casevm5a_case60t'#'../check_code/case11a'
#ddir = './'
#ddeg=0.25
#nlon=int(360/ddeg)+1
#nlat=int(180/ddeg)+1
#lon = np.arange(0,360.+ddeg,ddeg)
#lat = 90.-np.arange(0,180.+ddeg,ddeg)

ddeg=1.
nlon=int(360/ddeg)
nlat=int(180/ddeg)
#pdb.set_trace()
#lon = np.arange(0,360.,ddeg)
#lat = 90.-np.arange(0.5,180.5,ddeg)

#in CitcomSVE the lon is 0 ~ 360; lat is -90 ~ 90
lon = np.arange(0,361.,ddeg)
lat = np.arange(0,181.,ddeg)-90.

clm = np.zeros((lmax+1,lmax+1))
slm = np.zeros((lmax+1,lmax+1))
uclm = np.zeros((lmax+1,lmax+1))
uslm = np.zeros((lmax+1,lmax+1))
hclm = np.zeros((lmax+1,lmax+1))
hslm = np.zeros((lmax+1,lmax+1))
plm = plm_mohlenkamp(lmax,np.sin(lat*np.pi/180.0))

#caseid = 'PW.ice6g.1'+casetime
caseid = 'PW.ice6g.'+casetime
infile=caseid
outfile1='grate_grid_'+caseid+'.nc'
outfile2='urate_grid_'+caseid+'.nc'
outfile3='erate_grid_'+caseid+'.nc'
outfile4='nrate_grid_'+caseid+'.nc'
outfile5='rslrate_grid_'+caseid+'.nc'

data = np.loadtxt(path.join(ddir,infile),dtype=np.float)
# data format is (l, m, clm, slm, uclm, uslm, hclm, hslm), i.e. degree, order, geoid stokes coeff, uplift Stokes, and horizontal Stokes

i=0
for l in range(lmax+1):
        #print(l)
	for m in range(l+1):
		#if (l!=data[i,0]) or (m!=data[i,1]):
		#	print "something is wrong!", l, m, data[i,0], data[i,1]
		#	stop
                if m==0:
                        coef = 1./np.sqrt(4*M_PI)
                else:
                        coef = 1./(np.sqrt(2*M_PI)*(-1)**m)

		clm[l,m] = data[i,2]*coef
		slm[l,m] = -1*data[i,3]*coef
		uclm[l,m] = data[i,4]*coef
		uslm[l,m] = -1*data[i,5]*coef
		hclm[l,m] = data[i,6]*coef
		hslm[l,m] = -1*data[i,7]*coef
		i = i+1

grate = calc_sig(clm,slm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
urate = calc_sig(uclm,uslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
erate,nrate = calc_hori(hclm,hslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)
erate = erate*fact
nrate = nrate*fact

ncdf_write(grate.T,lon,lat,0,FILENAME=path.join(ddir,outfile1))
ncdf_write(urate.T,lon,lat,0,FILENAME=path.join(ddir,outfile2))
ncdf_write((grate-urate).T,lon,lat,0,FILENAME=path.join(ddir,outfile5))
ncdf_write(erate.T,lon,lat,0,FILENAME=path.join(ddir,outfile3))
ncdf_write(nrate.T,lon,lat,0,FILENAME=path.join(ddir,outfile4))

# save text files:
if 0:
	otxt1='grate_stokes_'+caseid+'.txt'
	np.savetxt(otxt1,data[:,:4],fmt='   %3d   %3d %14.5e  %14.5e  ')

if 0:
        otxt2='urate_grid_'+caseid+'.txt'
	f = open(otxt2,'w')
# lon:0...360, lat:90...(-90), Uplift, East, North
	for p in range(len(lat)):
		for q in range(len(lon)):
			f.write('%7.2f   %7.2f   %14.5f  %14.5f  %14.5f\n'%(lon[q],lat[p],urate[q,p],erate[q,p],nrate[q,p]))
	f.close()

#generate the 1 by 1 degree grid geoid/uplift and use the scaling factor in CitcomSVE to output the dementionless value
#scaling factor used in CitcomSVE
fac_pi = np.sqrt(4*np.pi)
fac_u = 6370000
fac_u = 1./fac_u

if 1:
        #otxt3='u_grid_'+str(case_time)+'.txt'
	otxt3='CASE'+caseid0+'.tps_grid.'+casetime
        f = open(path.join(ddir,otxt3),'w')
# lon:0...360, lat:90...(-90), Uplift, East, North
        for p in range(len(lat)):
                for q in range(len(lon)):
                        f.write('%7.2e   %7.2e   %14.10e\n'%(lon[q],lat[p],urate[q,p]*fac_u))
        f.close()

#fac_g = 4*np.pi*6.67*10**(-11)*4400*(6370000)**2/9.8
fac_g = 4*np.pi*6.67*10**(-11)*4604.4*(6370000)**2/9.8
fac_g = 1./fac_g

if 1:
        #otxt4='g_grid_'+str(case_time)+'.txt'
	otxt4='CASE'+caseid0+'.pttl_grid.'+casetime
        f = open(path.join(ddir,otxt4),'w')
# lon:0...360, lat:90...(-90), Uplift, East, North
        for p in range(len(lat)):
                for q in range(len(lon)):
                        f.write('%7.2e   %7.2e   %14.10e\n'%(lon[q],lat[p],grate[q,p]*fac_g))
        f.close()



#generate spherical harmonic coefficients which can directly be compared with those in CitcomSVE 
if 1:
	#otxt5='u_sh_'+str(case_time)+'.txt'
	otxt5='CASE'+caseid0+'.tps_sharm.'+casetime
	f = open(path.join(ddir, otxt5),'w')
	for l in range(lmax+1):
		for m in range(l+1):
                	if m==0:
                        	coef1 = 1.
                	else:
                        	coef1 = (-1)**m
			#f.write('%3d  %3d  %14.10e   %14.10e\n'%(l,m,uclm[l,m]*fac_u*fact,uslm[l,m]*fac_u*fact))
                        f.write('%3d  %3d  %14.10e   %14.10e\n'%(l,m,uclm[l,m]*fac_u*fact*fac_pi*coef1,uslm[l,m]*fac_u*fact*fac_pi*coef1))
	f.close()

if 1:
        #otxt6='g_sh_'+str(case_time)+'.txt'
	otxt6='CASE'+caseid0+'.pttl_sharm.'+casetime
        f = open(path.join(ddir, otxt6),'w')
        for l in range(lmax+1):
                for m in range(l+1):
                        if m==0:
                                coef1 = 1.
                        else:
                                coef1 = (-1)**m
                        #f.write('%3d  %3d  %14.10e   %14.10e\n'%(l,m,clm[l,m]*fac_g*fact,slm[l,m]*fac_g*fact))
                        f.write('%3d  %3d  %14.10e   %14.10e\n'%(l,m,clm[l,m]*fac_g*fact*fac_pi*coef1,slm[l,m]*fac_g*fact*fac_pi*coef1))
        f.close()
