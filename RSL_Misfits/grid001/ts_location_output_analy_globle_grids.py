#!/usr/bin/env python
# This code reads the Stokes coefficients of radial (UP) and horizontal (East and North) displacement in meters
# Generate the RSL history for a given location

from os import path
from plm_mohlenkamp import plm_mohlenkamp
from gen_clms import gen_clms
from calc_sig import calc_sig
from calc_hori import calc_hori
import numpy as np
from os import path
import re
import matplotlib.pyplot as plt
import pdb
from scipy import interpolate

lmax=100  #32  #100
M_PI = np.pi
fact=1.e-2 # convert cm to meters

nepochs = 50 #55 for ANU; 50 for ice6g
caseid = 'vm5a_compr' #'case_anu07' #'case_anu09'
ddir = './' #'../ice6g_compressible_ANU.dir/results_anu12' #'../ice6g_incompressible_present_day_ice.dir/results_case60t' #'../ice6g_compressible_ANU.dir/results_anu07' #'../ice6g_incompressible_ANU.dir/results_120p83b'
#of = open(path.join(ddir, 'Lambeck_rsl_predictions_'+caseid+'.txt'),'w')

#read eustatic sea level from compressible case
rsl_C = np.loadtxt(path.join(ddir,'RSL_c_vm5a_compr.txt'))
#rsl_C = np.loadtxt(path.join(ddir,'RSL_c_case60t.txt')) #'RSL_c_vm5a_compr.txt'))

#generate the global grids
lats = np.arange(90.,-91.,-1.0) #np.arange(-90.0,91.0,1.0)
lons = np.arange(0.,361.,1.0)

rsl = np.zeros([len(lats)*len(lons),nepochs])*np.nan
uplift = np.zeros([len(lats)*len(lons),nepochs])*np.nan
geoid = np.zeros([len(lats)*len(lons),nepochs])*np.nan


for count in range(nepochs):
	print(count)

	plm = plm_mohlenkamp(lmax,np.sin(lats*np.pi/180.0))
	gclm = np.zeros((lmax+1,lmax+1))
	gslm = np.zeros((lmax+1,lmax+1))
	uclm = np.zeros((lmax+1,lmax+1))
	uslm = np.zeros((lmax+1,lmax+1))

	infile='PW.ice6g.'+str(101+count)
	data = np.loadtxt(path.join(ddir,infile),dtype=np.float)

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


			gclm[l,m] = data[i,2]*coef
			gslm[l,m] = -data[i,3]*coef
			uclm[l,m] = data[i,4]*coef
			uslm[l,m] = -data[i,5]*coef
			i = i+1

	grate = calc_sig(gclm,gslm,lons,lats,LMIN=0,LMAX=lmax,PLM=plm)*fact
	urate = calc_sig(uclm,uslm,lons,lats,LMIN=0,LMAX=lmax,PLM=plm)*fact
	rsl[:,count] = grate.T.reshape([361*181]) - urate.T.reshape([361*181]) + rsl_C[count]
	geoid[:,count] = grate.T.reshape([361*181])
	uplift[:,count] = urate.T.reshape([361*181])

	
#write to file
#plt.figure()
#plt.plot(rsl[0,:])
#plt.show()
#np.savetxt(path.join(ddir,'RSL_Lambeck_'+caseid+'.txt'),rsl)
np.savetxt(path.join(ddir,'RSL_ICE6G_'+caseid+'.txt'),rsl)
np.savetxt(path.join(ddir,'GEOID_ICE6G_'+caseid+'.txt'),geoid)
np.savetxt(path.join(ddir,'UPLIFT_ICE6G_'+caseid+'.txt'),uplift)

