#!/usr/bin/env python
# This code is to read Lambeck's far field RSL dataset

from os import path
import numpy as np
import matplotlib.pyplot as plt
import pandas
import pdb

ddir = "../"
#fn = "Lambeck_farfield_RSL/Lambeck_etal_2014_PNAS_supp_edit.xlsx"  # Far field
#data = pandas.read_excel(path.join(ddir, fn))
#sitename = data['Site'].values
#RSL = data['RSL'].values
#RSLerr = data['Sig.RSL'].values
#age = data['Age'].values
#lat = data['Latitude'].values
#lon = data['Longitude'].values
#sample = data['Sample'].values


fn = "ANU_ice_model/QSR_SI/rsl_observations/rsl data set (after Dyke 2004).xlsx" # North America
data = pandas.read_excel(path.join(ddir, fn))
sitename = data['site name'].values
RSL = data['observed rsl'].values
RSLerr = data['2xsigma rsl'].values
age = data['time ka BP'].values
lat = data['latitude N'].values
lon = data['longitude E'].values
#sample = data['Sample'].values

# Extract site name and remove the duplicates
sitename_list = []
for i in sitename:
        if i not in sitename_list:
                sitename_list.append(i)

lon_ave_list = []#np.zeros(len(UTKEY_list))
lat_ave_list = []#np.zeros(len(UTKEY_list))
sitename_ave_list = []  # save the sitenames while selected critira applied to the dataset, e.g. num. point > 5

count = 0
for sitename_i in sitename_list:
        print(sitename_i)
        #pdb.set_trace()
        # Set random integer
        #random_integer = np.random.randint(0, len(UTKEY_list))
        #UTKEY_rand = UTKEY_list[random_integer]

        #idx = (UTKEY == UTKEY_rand)
        idx = (sitename == sitename_i)

        # if the number of data is less than 5, then ignore the site
        site_rsl = RSL[idx]

        if len(site_rsl) <= 5:
                continue

        site_rsl_err = RSLerr[idx]
        site_age = age[idx]
        site_lon = lon[idx]
        site_lat = lat[idx]

        # Detect and delete outliers
        mean_site_rsl  = np.mean(site_rsl)
        std_site_rsl = np.std(site_rsl)
        idx_del_outlier = ((site_rsl >= mean_site_rsl - 3*std_site_rsl) & (site_rsl <= mean_site_rsl + 3*std_site_rsl))

        site_rsl1 = site_rsl[idx_del_outlier]
        site_rsl_err1 = site_rsl_err[idx_del_outlier]
        site_age1 = site_age[idx_del_outlier]
        site_lon1 = site_lon[idx_del_outlier]
        site_lat1 = site_lat[idx_del_outlier]

        # Calculate the average of lon, lat for each utkey
        lon_ave_list.append(np.mean(site_lon1))
        lat_ave_list.append(np.mean(site_lat1))
        sitename_ave_list.append(sitename_i)

        plt.figure();
        plt.plot(-site_age,site_rsl,'.r')
        plt.errorbar(-site_age, site_rsl, yerr=site_rsl_err, fmt='.r', marker='', mfc='r', mec='r', ecolor='r')

        plt.plot(-site_age1,site_rsl1,'.k')
        plt.errorbar(-site_age1, site_rsl1, yerr=site_rsl_err1, fmt='.k', marker='', mfc='k', mec='k', ecolor='k')

        plt.xlabel('Time(kyr)')
        plt.ylabel('RSL(meters)')
        plt.title('Lambeck Far Field RSL at Site: {}'.format(sitename_i))
        #plt.show()
        plt.close()

        count += 1

print(count)

#pdb.set_trace()
# Write the longitude and latitude to file and save
lon_ave_arr = np.array(lon_ave_list)
lat_ave_arr = np.array(lat_ave_list)

# Remove the data with averaged latitude > 90.0 (wrong data)
idx_lat_ave_arr = (lat_ave_arr >= -90.) & (lat_ave_arr <= 90.)
lat_ave_arr = lat_ave_arr[idx_lat_ave_arr]
lon_ave_arr = lon_ave_arr[idx_lat_ave_arr]

# The longitude range in Lambeck dataset is 0 to 360.
# Remove the data with averaged longitude > 360.0 (wrong data)
idx_lon_ave_arr = (lon_ave_arr >= 0.0) & (lon_ave_arr <= 360.)
lat_ave_arr = lat_ave_arr[idx_lon_ave_arr]
lon_ave_arr = lon_ave_arr[idx_lon_ave_arr]

#np.savetxt('RSL_site_location_selected_Lambeck_FarField.txt', np.array([lon_ave_arr, lat_ave_arr]).T)
#np.savetxt('Sitename_Lambeck_FarField_selected.txt', sitename_ave_list, fmt='%s')

np.savetxt('RSL_site_location_selected_Lambeck_NA.txt', np.array([lon_ave_arr, lat_ave_arr]).T)
np.savetxt('Sitename_Lambeck_NA_selected.txt', sitename_ave_list, fmt='%s')


