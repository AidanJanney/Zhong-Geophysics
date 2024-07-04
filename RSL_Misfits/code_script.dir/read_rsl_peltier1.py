#!/usr/bin/env python
#this code is to read petlier's rsl excel table 
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas
import matplotlib
from scipy import interpolate


ddir = '../'
fn = 'rsl_dbase_peltier_all.xlsx'
 
data = pandas.read_excel(path.join(ddir,fn),sheet_name = 'dbase_table')
UTKEY = data['UTKEY'].values
TYP = data['TYP'].values
RSL = data['RSL'].values
RSLMIN = data['RSLMIN'].values
RSLMAX = data['RSLMAX'].values
AGE = data['AGE'].values
AGEMIN = data['AGEMIN'].values
AGEMAX = data['AGEMAX'].values
LAT = data['LAT'].values
LON = data['LON'].values


#range
LON1,LON2 = -110.0,-60.0 
LAT1,LAT2 = 50.0,68.0

idx_lon_r = (LON >= LON1) & (LON <= LON2)
idx_lat_r = (LAT >= LAT1) & (LAT <= LAT2)
idx_r = idx_lon_r & idx_lat_r

#crop the data in the range
UTKEY_r = UTKEY[idx_r]
TYP_r = TYP[idx_r]
RSL_r = RSL[idx_r]
RSLMIN_r = RSLMIN[idx_r]
RSLMAX_r = RSLMAX[idx_r]
AGE_r = AGE[idx_r]
AGEMIN_r = AGEMIN[idx_r]
AGEMAX_r = AGEMAX[idx_r]
LAT_r = LAT[idx_r]
LON_r = LON[idx_r]

#read data in the range by UTKEY, if the average location of the nearby UTKEY is less than 1 degree then combine
#extract UTKEY and remove the duplicates
UTKEY_list = []
for i in UTKEY_r:
        if i not in UTKEY_list:
                UTKEY_list.append(i)

LON_ave_arr = np.zeros(len(UTKEY_list))
LAT_ave_arr = np.zeros(len(UTKEY_list))
k=0
#pdb.set_trace()
#calculate the average location for each UTKEY_list
for i in UTKEY_list:
	#read the current UTKEY
	idx_i = (UTKEY_r == i)
	LON_i = LON_r[idx_i]
	LAT_i = LAT_r[idx_i]
	LON_i_ave = np.mean(LON_i)
	LAT_i_ave = np.mean(LAT_i)
	LON_ave_arr[k] = LON_i_ave
	LAT_ave_arr[k] = LAT_i_ave
	k += 1
#UTKEY_list = np.array(UTKEY_list)

#group data in the UTKEY_list by location
lim_lon = 0.5 #0.5  #degree
lim_lat = 0.5 #0.5  #degree
UTKEY_g_id = []

#inicialize array
#UTKEY_list_new = UTKEY_list
LON_ave_list = LON_ave_arr.tolist()
LAT_ave_list = LAT_ave_arr.tolist()

#pdb.set_trace()
for i in np.zeros(len(UTKEY_list)):
	i = int(i)
	LON_ave_arr = np.array(LON_ave_list)
	LAT_ave_arr = np.array(LAT_ave_list)
	#UTKEY_list = UTKEY_list_new

	#LON_ave_list = LON_ave_arr.tolist()
	
	idx_lon_g = (LON_ave_arr >= LON_ave_arr[i] - lim_lon) & (LON_ave_arr <= LON_ave_arr[i] + lim_lon)	
	idx_lat_g = (LAT_ave_arr >= LAT_ave_arr[i] - lim_lat) & (LAT_ave_arr <= LAT_ave_arr[i] + lim_lat)
	idx_g = idx_lon_g & idx_lat_g

	UTKEY_g = np.array(UTKEY_list)[idx_g]
	#test
        #print(UTKEY_g)
        #print(UTKEY_list)
	UTKEY_g_id.append(UTKEY_g)
	UTKEY_list = [ele for ele in UTKEY_list if ele not in UTKEY_g]

	#test
	#print(i)
	#print(UTKEY_g)
	#print(UTKEY_list)

	LON_ave_list = LON_ave_arr.tolist()
	LON_ave_list = [ele for ele in LON_ave_list if ele not in LON_ave_arr[idx_g]]

	LAT_ave_list = LAT_ave_arr.tolist()
	LAT_ave_list = [ele for ele in LAT_ave_list if ele not in LAT_ave_arr[idx_g]]

	if UTKEY_list == []:
		break;
#pdb.set_trace()
#plot 

for UTKEY_i in UTKEY_g_id:
	plt.figure();
	print(UTKEY_i)
	for UTKEY_ii in UTKEY_i:
		idx_UTKEY_ii = (UTKEY == UTKEY_ii)
		
		site_rsl = RSL[idx_UTKEY_ii]
		site_rsl_min = RSLMIN[idx_UTKEY_ii]
		site_rsl_max = RSLMAX[idx_UTKEY_ii]
		site_age = AGE[idx_UTKEY_ii]/1000.
		site_age_min = AGEMIN[idx_UTKEY_ii]/1000.
		site_age_max = AGEMAX[idx_UTKEY_ii]/1000.

		plt.plot(-site_age,site_rsl,'.k')
		plt.errorbar(-site_age,site_rsl,xerr=np.array([site_age_max-site_age,site_age-site_age_min]),yerr=np.array([site_rsl-site_rsl_min,site_rsl_max-site_rsl]),fmt='.k', marker='', mfc='k',mec='k', ecolor='k')

	plt.show()


