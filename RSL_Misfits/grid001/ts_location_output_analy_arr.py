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
import matplotlib.pyplot as plt
import pdb
from scipy import interpolate

lmax=100  #32  #100
M_PI = np.pi
fact=1.e-2 # convert cm to meters

ddir_of = './' #'../check_code/casevm5a_case60t_50epoch/ts_sites_output.dirr' #'../check_code/case_incompr_Lambeck/ts_sites_output.dirr'#'../ice6g_incompressible_ANU.dir/results_anu01a/ts_sites_output.dirr' #'../ice6g_incompressible_ANU.dir/results_50p12/ts_sites_output.dirr'#'../check_code/casevm5a_case60t_50epoch/ts_sites_output.dirr' #'../ice6g_incompressible_ANU.dir/results_50p12/ts_sites_output.dirr'#'../check_code/casevm5a_case60t_50epoch/ts_sites_output.dirr' #'../ice6g_incompressible_ANU.dir/results_50p12/ts_sites_output.dirr'
nepochs = 50 #245 #55 for ANU #93#50 #55 for W12#114 for Lambeck #122 for ice6g #50 for ice6g
#casetag = 'case_compr_' #'case_compr_' #'CASE'
caseid = 'casevm5a_compr'#'case_anu01a'#'case60t' #'case_incompr_vm5a_shear_ocean120_Glenn'#'case_incompr_wei' #'case_incompr_W12_b'#'case_incompr_Lambeck2'#'case_incompr_ICE7g_vm7'#'case_incompr_W12_a' #'case_incompr_Lambeck_global_high' #'case_incompr_Lambeck' #'case_incompr_Lambeck_global_high' #'case_incompr_Lambeck1'#'case_incompr_Lambeck' #'case_compr_vm5a' #'01' #'vm5a' #'30fac'    #case11: 11c vs 11; case12: 12c vs 12a;  case13: 13a vs 13

#sitename = 'Barbados'
#lon = np.array([360-59.53, 360-59.53])
#lat = np.array([13.07, 13.07])

#sitename = 'Churchill'
#lon = np.array([265.6, 265.6])
#lat = np.array([58.7, 58.7])

#sitename = 'Richmond_Gulf'
#lon = np.array([283.6, 283.6])
#lat = np.array([56.6, 56.6])

#sitename = 'Boston'
#lon = np.array([360-71.1,360-71.1])
#lat = np.array([42.4, 42.4])

#sitename = 'Mississippi'
#lon = np.array([269.8, 269.8])
#lat = np.array([30.0, 30.0])

#sitename = 'Sotra_Norway'
#lon = np.array([5.1, 5.1])
#lat = np.array([60.3, 60.3])

#sitename = 'Pt_Caen'
#lon = np.array([360-115.8, 360-115.8])
#lat = np.array([69.3, 69.3])

#sitename = 'Murray_Bay'
#lon = np.array([360-93.7, 360-93.7])
#lat = np.array([71.7, 71.7])

#site in Maine
sitename_arr1 = ['Site 1','Site 2','Site 3','Site 4','Site 5','Site 6','Site 7a','Site 7b','Site 8a','Site 8b','Site 9','Site 10']
loc_arr1 = np.array([[360.-64.203,45.609],[360.-63.383,44.639],[360.-60.025,43.945],[360.-65.618,43.728],[360.-65.903,44.530],[360.-66.606,45.104],[360.-68.836,44.968],[360.-67.809,44.696],[360.-69.917,44.474],[360.-70.074,43.806],[360.-70.801,42.702],[360.-70.527,41.579]])

sitename_arr2 = ['F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12']
loc_arr2 = np.array([[29.4,70.1],[19.0,69.8],[10.9,63.6],[5.1,60.3],[9.4,58.9],[10.8,59.8],[12.0,57.4],[17.7,59.5],[19.9,64.0],[23.3,60.1],[24.9,59.4],[10.7,57.7]])

sitename_arr3 = ['N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18']
loc_arr3 = np.array([[304.2,51.5],[295.7,59.8],[290.0,59.4],[285.6,62.1],[282.1,53.6],[265.6,58.7],[244.9,67.9],[232.0,52.0],[236.8,49.6],[263.2,28.2],[269.8,30.0],[279.8,25.4],[279.6,32.7],[286.1,41.2],[289.5,41.6],[289.9,43.6],[295.7,45.8],[290.5,47.8]])

sitename_arr4 = ['Tahiti','Kiritimati','Christchurch','Port Pirie','Redcliff Belperio','Port Gawler','Fisherman Bay','Wood Point','Abrolhos','Grub Reef','Pioneer Bay','Magnetic Island','Yule Point','Karumba','Huon Peninsula','Semakau','Geylang','Phuket','Ca_Na','Maldives Maalhosmadulu','Maldives Rasdhoo','Reunion','Seychelles','Mayotte','Natal','Barbados']
loc_arr4 = np.array([[210.45,-17.45],[202.52,1.98],[172.67,-43.5],[138.01,-33.16],[137.87,-32.69],[138.46,-34.64],[137.88,-33.43],[137.86,-33.33],[113.83,-28.68],[147.43,-18.63],[146.48,-18.60],[146.87,-19.15],[145.52,-16.57],[140.83,-17.42],[147.62,-6.11],[103.77,1.21],[103.87,1.31],[98.41,7.75],[108.84,11.33],[73.03,5.27],[72.98,4.30],[55.23,-21.08],[55.52,-4.68],[45.27,-12.80],[32.40,-28.47],[300.45,13.04]])

#sitename_arr = ['N5','F9','Barbados']
#loc_arr = np.array([[282.1,53.6],[19.9,64.0],[300.45,13.04]])

#sitename_arr = sitename_arr1 + sitename_arr2 + sitename_arr3 +sitename_arr4
sitename_arr5 = ['W Ungava Bay','Deception Bay','Kugluktuk','Riviere du Loup','Goulburn Lake','S Massachusetts','New York','Ft George','Frosta Norway','Ski Moraine Norway','Sandsjobacka Sweden','Kragero Norway','Sotra Norway','Borgan Norway','Vestfonna','Varangerfjord Norway']#,'Hinnerjoki Finland']#'Vestfonna','Spitsbergen']
loc_arr5 = np.array([[290.0,59.4],[285.6,62.1],[244.9,67.9],[290.5,47.8],[252.0,67.0],[289.5,41.6],[286.1,41.2],[282.1,53.6],[10.9,63.6],[10.8,59.8],[12.0,57.4],[9.4,58.9],[5.1,60.3],[11.0,65.0],[19.0,80.0],[29.4,70.1]])#[22.0,61.0]]) #,[19.0,80.0],[14.0,78.0]

sitename_arr6 = ['Barbados','Huon Peninsula','Mayotte']
loc_arr6 = np.array([[300.45,13.04],[147.62,-6.11],[45.27,-12.80]])

sitename_arr = sitename_arr2 + sitename_arr3 + sitename_arr4

#loc_arr = np.append(loc_arr1,loc_arr2,0)
loc_arr = np.append(loc_arr2,loc_arr3,0)
loc_arr = np.append(loc_arr,loc_arr4,0)
#loc_arr = np.append(loc_arr5,loc_arr6,0)

for i_site in range(len(sitename_arr)):
#for i_site in [0,1]:
	#i_site = 0#25
	print(i_site)
	sitename = sitename_arr[i_site]
	lon = np.array([loc_arr[i_site,0],loc_arr[i_site,0]])
	lat = np.array([loc_arr[i_site,1],loc_arr[i_site,1]])

	#sitename = sitename_arr
	#lon = loc_arr[:,0]
	#lat = loc_arr[:,1]

	plm = plm_mohlenkamp(lmax,np.sin(lat*np.pi/180.0))
	gclm = np.zeros((lmax+1,lmax+1))
	gslm = np.zeros((lmax+1,lmax+1))
	uclm = np.zeros((lmax+1,lmax+1))
	uslm = np.zeros((lmax+1,lmax+1))
	hclm = np.zeros((lmax+1,lmax+1))
	hslm = np.zeros((lmax+1,lmax+1))

	rsl = np.zeros(nepochs)*np.nan
	geoid = np.zeros(nepochs)*np.nan
	uplift = np.zeros(nepochs)*np.nan
	east = np.zeros(nepochs)*np.nan
	north = np.zeros(nepochs)*np.nan

	#output file
	of_rsl = sitename+'_rsl_'+caseid+'.txt'
	of_geoid = sitename+'_geoid_'+caseid+'.txt'
	of_uplift = sitename+'_uplift_'+caseid+'.txt'
	of_east = sitename+'_east_'+caseid+'.txt'
	of_north = sitename+'_north_'+caseid+'.txt'

	#read eustatic sea level
	#rsl_C = np.loadtxt('RSL_c.txt')

	#read eustatic sea level from compressible case
	ddir= '.' #'../check_code/casevm5a_case60t_50epoch' #'../check_code/case_incompr_Lambeck'#'../ice6g_incompressible_ANU.dir/results_anu01a' #'../ice6g_incompressible_ANU.dir/results_50p12'#'../check_code/casevm5a_case60t_50epoch' #'../ice6g_incompressible_ANU.dir/results_50p12' #'../check_code/casevm5a_case60t' #'../ice6g_incompressible_arbitr_time.dir/' #'../check_code/casevm5a_case60t'#'../ice6g_incompressible_ANU.dir/results_120p83b' #'../ice6g_compressible_ANU.dir/results_anu07' #'../ice6g_incompressible_ANU.dir/results_120p83_ice6g0' #'../ice6g_incompressible_ice_Glenn.dir/results_ocean_ice6g120p55_vm5a_shear.dir'#'../check_code/case_incompr_wei'#'../check_code/case_incompr_W12_b'#'../check_code/case_incompr_Lambeck2'#'../check_code/case_incompr_ICE7g_vm7/'#'../check_code/case_incompr_W12_a/'#'../check_code/case_incompr_Lambeck_global_high' #'../check_code/case_incompr_Lambeck_global_high' #'../check_code/case_incompr_Lambeck1'  #case_compr_Lambeck_122'  #case_incompr_Lambeck/' #'../ice6g_compressible.dir/'
	
	#rsl_C = np.loadtxt(path.join(ddir,'RSL_c_Glenn.txt'))  #ANU
	#rsl_C = np.loadtxt('./RSL_c_50epoch.txt')  #ICE6G
	rsl_C = np.loadtxt('./RSL_c_vm5a_compr.txt')  #ICE6G+vm5a+compr
	
#pdb.set_trace()
	#interpolate RSL_c of ice6g 122 to Lambeck ice epoch, if the latter one is calculated,then comment the following lines
	#rsl_time1 = np.loadtxt(path.join(ddir,'time.reverse'))
	#rsl_time2 = np.loadtxt('./years.reverse')
	#rsl_f = interpolate.interp1d(rsl_time2,rsl_C,fill_value='extrapolate')
	#rsl_C1 = rsl_f(rsl_time1)
	#rsl_C = rsl_C1

	for count in range(nepochs):
		#ddir='../ice6g_compressible.dir/'
		infile='PW.ice6g.'+str(101+count)
		data = np.loadtxt(path.join(ddir,infile),dtype=np.float)
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

				#exclude the 2,1 terms
				#if l == 2 and m == 1:
				#	coef = 0.0

				gclm[l,m] = data[i,2]*coef
				gslm[l,m] = -data[i,3]*coef
				uclm[l,m] = data[i,4]*coef
				uslm[l,m] = -data[i,5]*coef
                		hclm[l,m] = data[i,6]*coef
                		hslm[l,m] = -data[i,7]*coef
				i = i+1

		grate = calc_sig(gclm,gslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
		urate = calc_sig(uclm,uslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)*fact
		erate,nrate = calc_hori(hclm,hslm,lon,lat,LMIN=0,LMAX=lmax,PLM=plm)
		erate = erate*fact
		nrate = nrate*fact

		rsl[count] = grate[0,0] - urate[0,0]
		geoid[count] = grate[0,0]
		uplift[count] = urate[0,0]
		east[count] = erate[0,0]
		north[count] = nrate[0,0]

	print grate.shape

	#read years.reverse
	time_year = np.loadtxt(path.join(ddir,'years.reverse'))
	#time_year = np.loadtxt('time.reverse')
	#time_year = np.loadtxt('years.reverse1')  #50 epoch
	time = -1*time_year

	#pdb.set_trace()
	f1 = open(path.join(ddir_of, of_rsl),'w')
	f2 = open(path.join(ddir_of, of_uplift),'w')
	f3 = open(path.join(ddir_of, of_geoid),'w')
	f4 = open(path.join(ddir_of, of_east),'w')
	f5 = open(path.join(ddir_of, of_north),'w')


	rsl = rsl + rsl_C

	for i in range(nepochs):
		f1.write('%14.10e   %14.10e\n'%(time[i],rsl[i]))
		f2.write('%14.10e   %14.10e\n'%(time[i],uplift[i]))
		f3.write('%14.10e   %14.10e\n'%(time[i],geoid[i]))
		f4.write('%14.10e   %14.10e\n'%(time[i],east[i]))
		f5.write('%14.10e   %14.10e\n'%(time[i],north[i]))

	f1.close()
	f2.close()
	f3.close()
