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
fn = 'RSL_measurements/rsl_Peltier_2015paper.xls'
 
data = pandas.read_excel(path.join(ddir,fn),sheet_name = 'sheet 1')
UTKEY = data['UTKEY'].values
#TYP = data['TYP'].values
RSL = data['RSL'].values
RSLMIN = data['RSLMIN'].values
RSLMAX = data['RSLMAX'].values
AGE = data['AGE'].values
AGEMIN = data['AGEMIN'].values
AGEMAX = data['AGEMAX'].values
LAT = data['LAT'].values
LON = data['LON'].values

#read data in the range by UTKEY, if the average location of the nearby UTKEY is less than 1 degree then combine
#extract UTKEY and remove the duplicates
UTKEY_list = []
for i in UTKEY:
        if i not in UTKEY_list:
                UTKEY_list.append(i)

LON_ave_list = []#np.zeros(len(UTKEY_list))
LAT_ave_list = []#np.zeros(len(UTKEY_list))
UTKEY_ave_list = []

count = 0
for UTKEY_i in UTKEY_list:
	print(UTKEY_i)
	#pdb.set_trace()
	# Set random integer
	#random_integer = np.random.randint(0, len(UTKEY_list))
	#UTKEY_rand = UTKEY_list[random_integer]
	
	#idx = (UTKEY == UTKEY_rand)
	idx = (UTKEY == UTKEY_i)

	# if the number of data is less than 5, then ignore the site
	site_rsl = RSL[idx]

	if len(site_rsl) <= 5:
		continue

	site_rsl_min = RSLMIN[idx]
	site_rsl_max = RSLMAX[idx]

	# ONLY applied Peltier's global rsl database
	#site_age = AGE[idx]/1000.
	#site_age_min = AGEMIN[idx]/1000.
	#site_age_max = AGEMAX[idx]/1000.

	# Applied on Peltier's digitized database and Lambeck's NA & farfield rsl database
        site_age = AGE[idx]
        site_age_min = AGEMIN[idx]
        site_age_max = AGEMAX[idx]
	site_lon = LON[idx]
	site_lat = LAT[idx]

        # Detect and delete outliers
        mean_site_rsl  = np.mean(site_rsl)
        std_site_rsl = np.std(site_rsl)
	idx_del_outlier = ((site_rsl >= mean_site_rsl - 3*std_site_rsl) & (site_rsl <= mean_site_rsl + 3*std_site_rsl))

	site_rsl1 = site_rsl[idx_del_outlier]
	site_rsl_min1 = site_rsl_min[idx_del_outlier]
	site_rsl_max1 = site_rsl_max[idx_del_outlier]
	site_age1 = site_age[idx_del_outlier]
	site_age_min1 = site_age_min[idx_del_outlier]
	site_age_max1 = site_age_max[idx_del_outlier]
	site_lon1 = site_lon[idx_del_outlier]
	site_lat1 = site_lat[idx_del_outlier]	

	# Calculate the average of lon, lat for each utkey
	LON_ave_list.append(np.mean(site_lon1))
	LAT_ave_list.append(np.mean(site_lat1))  
	UTKEY_ave_list.append(UTKEY_i)

	plt.figure();
	plt.plot(-site_age,site_rsl,'.r')
	plt.errorbar(-site_age,site_rsl,xerr=np.array([site_age_max-site_age,site_age-site_age_min]),yerr=np.array([site_rsl-site_rsl_min,site_rsl_max-site_rsl]),fmt='.r', marker='', mfc='r',mec='r', ecolor='r')

	plt.plot(-site_age1,site_rsl1,'.k')
	plt.errorbar(-site_age1,site_rsl1,xerr=np.array([site_age_max1-site_age1,site_age1-site_age_min1]),yerr=np.array([site_rsl1-site_rsl_min1,site_rsl_max1-site_rsl1]),fmt='.k', marker='', mfc='k',mec='k', ecolor='k')
	
	plt.xlabel('Time(kyr)')
	plt.ylabel('RSL(meters)')
	plt.title('Peltier RSL at Site: {}'.format(UTKEY_i))
	#plt.show()
	plt.close()

	count += 1

print(count)

#pdb.set_trace()
# Write the longitude and latitude to file and save
Lon_ave_arr = np.array(LON_ave_list)
Lat_ave_arr = np.array(LAT_ave_list)

# Remove the data with averaged latitude > 90.0 (wrong data)
#idx_Lat_ave_arr = (Lat_ave_arr >= -90.) & (Lat_ave_arr <= 90.)
#Lat_ave_arr = Lat_ave_arr[idx_Lat_ave_arr]
#Lon_ave_arr = Lon_ave_arr[idx_Lat_ave_arr]

# Remove the data with averaged longitude > 180.0 (wrong data)
#idx_Lon_ave_arr = (Lon_ave_arr >= -180.) & (Lon_ave_arr <= 180.)
#Lat_ave_arr = Lat_ave_arr[idx_Lon_ave_arr]
#Lon_ave_arr = Lon_ave_arr[idx_Lon_ave_arr]

np.savetxt('RSL_site_location_digitized_selected.txt', np.array([Lon_ave_arr, Lat_ave_arr]).T)
np.savetxt('UTKEY_pelier_digitized_selected.txt', UTKEY_ave_list, fmt='%s')
