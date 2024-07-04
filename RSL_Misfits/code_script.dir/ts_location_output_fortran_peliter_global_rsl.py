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
from scipy import interpolate
#import matplotlib.pyplot as plt
import pdb

#def calc_rsl_prediction(idGridFile):
#        print(idGridFile)
lmax=100  #32  #100
M_PI = np.pi
fact=1.e-2 # convert cm to meters

nepochs = 50 #55 for ANU; 50 for ice6g
caseid = 'grids_compr_peltier_rsl' #'grids_compr_lambeck_farfield' #'grids_compr_peltier_digitized'#'grids_compr_lambeck_NA' #'grids_compr_lambeck_farfield' #'grids_compr'
#ddir = '../grid' + idGridFile #'../ice6g_compressible_ANU.dir/results_anu12' #'../ice6g_incompressible_present_day_ice.dir/results_case60t' #'../ice6g_compressible_ANU.dir/results_anu07' #'../ice6g_incompressible_ANU.dir/results_120p83b'
#of = open(path.join(ddir, 'Lambeck_rsl_predictions_'+caseid+'.txt'),'w')
ddir = "./"
#read eustatic sea level from compressible case
rsl_C = np.loadtxt(path.join(ddir,'RSL_c_vm5a_compr.txt'))
#rsl_C = np.loadtxt(path.join(ddir,'RSL_c_case60t.txt')) #'RSL_c_vm5a_compr.txt'))

# Read longitude/latitude file
data_locations = np.loadtxt('/home/khuan/gridresearch_results_ICE6G_batch/RSL_site_peltier_location_selected.txt') 
lats = data_locations[:,1]
lons = data_locations[:,0]
lons[lons < 0] += 360.    #Peltier's longitude range is -180 to 180

#data_locations = np.loadtxt('/home/khuan/gridresearch_results_ANU_batch/RSL_site_location_selected_Lambeck_FarField.txt') 
#data_locations = np.loadtxt('/home/khuan/gridresearch_results_ANU_batch/RSL_site_location_selected_Lambeck_NA.txt') 
#data_locations = np.loadtxt('RSL_site_location_digitized_selected.txt') 
#lats = data_locations[:,1]
#lons = data_locations[:,0]  #Lambeck's longitude range is 0 to 360.

rsl = np.zeros([len(lats), nepochs])*np.nan
uplift = np.zeros([len(lats), nepochs])*np.nan
geoid = np.zeros([len(lats), nepochs])*np.nan


for count in range(nepochs):
	#print(count)

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
	#pdb.set_trace()
	grate = calc_sig(gclm,gslm,lons,lats,LMIN=0,LMAX=lmax,PLM=plm)*fact
	urate = calc_sig(uclm,uslm,lons,lats,LMIN=0,LMAX=lmax,PLM=plm)*fact
	rsl[:,count] = np.diagonal(grate) - np.diagonal(urate) + rsl_C[count]
	geoid[:,count] = np.diagonal(grate)
	uplift[:,count] = np.diagonal(urate)

	
#write to file
#plt.figure()
#plt.plot(rsl[0,:])
#plt.show()
#np.savetxt(path.join(ddir,'RSL_Lambeck_'+caseid+'.txt'),rsl)
np.savetxt(path.join(ddir,'RSL_ICE6G_'+caseid+'.txt'),rsl)
np.savetxt(path.join(ddir,'GEOID_ICE6G_'+caseid+'.txt'),geoid)
np.savetxt(path.join(ddir,'UPLIFT_ICE6G_'+caseid+'.txt'),uplift)

# Generate index for grided files
#idGridFileArr = np.arange(1, 43, 1).astype(str)

#for idx_grid in range(len(idGridFileArr)):
#        idGridFile = idGridFileArr[idx_grid].zfill(3)

# Call the function
#calc_rsl_prediction(idGridFile=idGridFile)
#        calc_rsl_prediction(idGridFile="_test")

