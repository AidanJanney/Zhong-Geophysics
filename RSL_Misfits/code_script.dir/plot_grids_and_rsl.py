#!/usr/bin/evn python
# This code is to plot all grids results against RSL data to check the potenial problems in the model setup

from os import path
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas

def plot_grids_per_site(sitename):
    """
    plot all the grids results against RSL data at a given site name
    """

    # Generate index for grided files
    idGridFileArr = np.arange(1, 43, 1).astype(str)

    # Find the index for the sitename in the selected sitename file
    ddir_sitename = '.'
    #fn_sitename = 'Sitename_Lambeck_FarField_selected.txt'
    fn_sitename = 'Sitename_Lambeck_NA_selected.txt'
    f_sitename = open(path.join(ddir_sitename, fn_sitename), 'r')
    selected_sitename = f_sitename.readlines()
    selected_sitename = [ele.strip() for ele in selected_sitename]
    idx_site = selected_sitename.index(sitename)

    plt.figure()
    for idx_case in range(len(idGridFileArr)):
	    ddirGridFile = './grid'+idGridFileArr[idx_case].zfill(3)
            #print(ddirGridFile)
	    #fnGridFile = 'RSL_ICE6G_grids_compr_lambeck_farfield.txt'
            fnGridFile = 'RSL_ICE6G_grids_compr_lambeck_NA.txt'
            rslModel = np.loadtxt(path.join(ddirGridFile, fnGridFile))
            rslModelAtSite = rslModel[idx_site, :]
            rslModelAtSite = rslModelAtSite - rslModelAtSite[-1]

            # Read time tag for RSL gridded model prediction
            ddirTime = '.'
	    fnTime = 'years.reverse'
	    rslModelTime = np.loadtxt(path.join(ddirTime, fnTime))

            TimeCutoff = 41.0 #kypb
            idx_LGM = (rslModelTime <= TimeCutoff)
            
            color = "grey"
            linewidth = "1."

            idx_minimum = [40, 41, 42]#[20, 30, 33, 24, 25, 40, 41, 42] # Far field: have local minimum rms at those grids
            idx_minimum = [23, 24, 25]#[20, 30, 33, 24, 25, 40, 41, 42] # Far field: have local minimum rms at those grids
            #if idx_case in idx_minimum:
            #    color = "blue"
            #    linewidth = "1.5"

            #idx_group = [26, 27, 1, 2, 3, 4, 5]
            #idx_group = [28, 29, 6, 7, 8, 9, 10]
            #idx_group = [30, 31, 11, 12, 13, 14, 15]
            #idx_group = [33, 32, 16, 17, 18, 19, 20]
            #idx_group = [34, 35, 21, 22, 23, 24, 25]
            idx_group = [36, 37, 38, 39, 40, 41, 42]

            if idGridFileArr[idx_case] == str(idx_group[0]):
            #if int(idGridFileArr[idx_case]) in [26, 27, 1, 2, 3, 4, 5]:
                color = "black"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[1]):
            #if int(idGridFileArr[idx_case]) in [28, 29, 6, 7, 8, 9, 10]:
                color = "brown"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[2]):
            #if int(idGridFileArr[idx_case]) in [30, 31, 11, 12, 13, 14, 15]:
                color = "blue"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[3]):
            #if int(idGridFileArr[idx_case]) in [33, 32, 16, 17, 18, 19, 20]:
                color = "cyan"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[4]):
            #if int(idGridFileArr[idx_case]) in [34, 35, 21, 22, 23, 24, 25]:
                color = "magenta"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[5]):
            #if int(idGridFileArr[idx_case]) in [36, 37, 38, 39, 40, 41, 42]:
                color = "yellow"
                linewidth = "1.5"

            if idGridFileArr[idx_case] == str(idx_group[6]):
                color = "green"
                linewidth = "1.5"

            plt.plot(-rslModelTime[idx_LGM], rslModelAtSite[idx_LGM], color=color, linewidth=linewidth)

    # Read in the RSL data at sitename
    # Far field
    #ddirRSL = '../'
    #fnRSL = 'Lambeck_farfield_RSL/Lambeck_etal_2014_PNAS_supp_edit.xlsx'
    #data = pandas.read_excel(path.join(ddirRSL, fnRSL))
    #SiteName = data['Site'].values
    #RSL = data['RSL'].values
    #RSLerr = data['Sig.RSL'].values
    #age = data['Age'].values
    #lat = data['Latitude'].values
    #lon = data['Longitude'].values
    #sample = data['Sample'].values

    # North America
    ddirRSL = '../'
    fnRSL = 'ANU_ice_model/QSR_SI/rsl_observations/rsl data set (after Dyke 2004).xlsx'
    data = pandas.read_excel(path.join(ddirRSL, fnRSL))
    SiteName = data['site name'].values
    RSL = data['observed rsl'].values
    RSLerr = data['2xsigma rsl'].values
    age = data['time ka BP'].values
    lat = data['latitude N'].values
    lon = data['longitude E'].values

    # Search rsl data for idx_site
    idx_rsl_site = (SiteName == sitename)

    rslDataAtSite = RSL[idx_rsl_site]
    rslerrAtSite = RSLerr[idx_rsl_site]
    ageAtSite = age[idx_rsl_site]

    # Plot RSL data together with gridded model predictions
    plt.errorbar(-ageAtSite, rslDataAtSite, yerr=rslerrAtSite, fmt='.r', marker='.', mfc='r', mec='r', ecolor='r')

    plt.xlim([-20, 0])
    #plt.xticks(np.arange(-20,0.3,2.5),['-20','','-15','','-10','','-5','','0'])
    plt.margins(0.05,0.05)
    plt.title(sitename, fontsize=12)
    plt.xlabel('Time (kybp)',fontsize=12)
    plt.ylabel('RSL (meters)',fontsize=12)

    plt.show()
    plt.close()

# Main function start here

# Read in the selected sitename
ddir_sitename = '.'
#fn_sitename = 'Sitename_Lambeck_FarField_selected.txt'
fn_sitename = 'Sitename_Lambeck_NA_selected.txt'
f_sitename = open(path.join(ddir_sitename, fn_sitename), 'r')
selected_sitename = f_sitename.readlines()
selected_sitename = [ele.strip() for ele in selected_sitename]

for sitename in selected_sitename:
    plot_grids_per_site(sitename)
