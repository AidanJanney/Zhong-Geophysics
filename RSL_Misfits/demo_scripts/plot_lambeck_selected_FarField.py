'''
Plot the RSL data for far-field sites from Lambeck
files:
        ../../RSL_OBS/RSL_dataset_Lambeck_farfield/Lambeck_etal_2014_PNAS_supp_edit.xlsx
'''

from os import path
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import pandas


ddirRSL = '../../RSL_OBS/RSL_dataset_Lambeck_farfield/'
fnRSL = 'Lambeck_etal_2014_PNAS_supp_edit.xlsx'

data = pandas.read_excel(path.join(ddirRSL, fnRSL))
sitename = data['Site'].values
RSL = data['RSL'].values
RSLerr = data['Sig.RSL'].values
age = data['Age'].values
lat = data['Latitude'].values
lon = data['Longitude'].values
sample = data['Sample'].values

# check site name
sitename_unique = np.unique(sitename) 
print('sitename_unique = ', sitename_unique)

Output_dir = './fig_Lambeck_FarField/'
if not path.exists(Output_dir):
    os.makedirs(Output_dir)

def plot_RSL_curve(age,RSL, RSLerr, title, fig_name):
    '''
    Plot RSL curve
    '''
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.errorbar(-1*age, RSL, yerr=RSLerr, fmt='o', color='k', ecolor='gray', elinewidth=3, capsize=0)
    ax.set_xlabel('Age (ka)')
    if(np.max(age) < 26):
        ax.set_xlim([-26,0])
    else:
        ax.set_xlim([-np.max(age),0])
    ax.set_ylabel('RSL (m)')
    ax.set_title(title)
    ax.grid(True)
    plt.savefig(path.join(Output_dir, fig_name))
    plt.close()

# Plot RSL curve for each site
for i in range(len(sitename_unique)):
    idx = np.where(sitename == sitename_unique[i])
    age_site = age[idx]
    RSL_site = RSL[idx]
    RSLerr_site = RSLerr[idx]
    lat_site = np.mean(lat[idx])
    lon_site = np.mean(lon[idx])
    sample_site = sample[idx]
    title = sitename_unique[i] + ' (lat, lon) = (' + str(lat_site) + ', ' + str(lon_site) + ')'
    fig_name = sitename_unique[i] + '.png'
    plot_RSL_curve(age_site,RSL_site, RSLerr_site, title, fig_name)
    print('plotting ', sitename_unique[i])