#!/usr/bin/evn python
# This code is to read RSL data from Lambeck's RSL database
# Read the gridded model predictions
# And plot the curve together to compare
# Then calculate the global averaged RMS

from os import path
import numpy as np
#import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import pandas

def calc_RMS_grid(idGridFile):
	''' 
        calcuale RMS for each grided file
        '''
	# Read RSL gridded model prediction
	#idGridFile = '001'
	ddirGridFile = '/home/khuan/ICE6G_results_0.1grids/grid'+idGridFile
	fnGridFile = 'RSL_ICE6G_grids_compr_lambeck_farfield.txt'
	#fnGridFile = 'RSL_ICE6G_grids_compr_lambeck_NA_ALL.txt'
	rslModel = np.loadtxt(path.join(ddirGridFile, fnGridFile))

	# Read time tag for RSL gridded model prediction
	ddirTime = '/home/khuan/ICE6G_results_0.1grids/grid001'
	fnTime = 'years.reverse'
	rslModelTime = np.loadtxt(path.join(ddirTime, fnTime))

	# Read RSL database from Lambeck based on selected sitename from file 'Sitename_Lambeck_FarField_selected.txt'
	# In the selected sitename file, the site which is less than 5 data points has been deleted.
	ddirRSL = '/home/khuan/rsl_database/'
	
        # Lamebek 2014 far field RSL database
	fnRSL = 'Lambeck_etal_2014_PNAS_supp_edit.xlsx'

	data = pandas.read_excel(path.join(ddirRSL, fnRSL))
	sitename = data['Site'].values
	RSL = data['RSL'].values
	RSLerr = data['Sig.RSL'].values
	age = data['Age'].values
	lat = data['Latitude'].values
	lon = data['Longitude'].values
	sample = data['Sample'].values

        # Lambeck 2017 North America RSL database
        #fnRSL = "rsl_dataset_lambeck_NA.xlsx" #'ANU_ice_model/QSR_SI/rsl_observations/rsl data set (after Dyke 2004).xlsx'
	#data = pandas.read_excel(path.join(ddirRSL, fnRSL))
	#sitename = data['site name'].values
	#RSL = data['observed rsl'].values
	#RSLerr = data['2xsigma rsl'].values
	#age = data['time ka BP'].values
	#lat = data['latitude N'].values
	#lon = data['longitude E'].values

	# Read selected sitename file
	ddir_sitename = '..'
	fn_sitename = 'Sitename_Lambeck_FarField_selected.txt'
	#fn_sitename = 'Sitename_Lambeck_NA_all.txt'
	#selected_sitename = np.loadtxt(path.join(ddir_sitename, fn_sitename), dtype=np.str) # loadtxt does not work here because the site name has space between words or uneven collums data
	f_sitename = open(path.join(ddir_sitename, fn_sitename), 'r')
	selected_sitename = f_sitename.readlines()


	err = np.zeros(len(selected_sitename))

	for idx_site in range(len(selected_sitename)):
	    #print(idx_site)
	    # Read model predictions
	    rslModelAtSite = rslModel[idx_site, :]
    	rslModelAtSite = rslModelAtSite - rslModelAtSite[-1]

	    # Search rsl data for idx_site
    	sitenameAtSite = selected_sitename[idx_site].strip() # delete /n in string
		idx_sitename = (sitenameAtSite == sitename)

		rslDataAtSite = RSL[idx_sitename]
		rslerrAtSite = RSLerr[idx_sitename]
		ageAtSite = age[idx_sitename]

	    # Make the plot for comparisons
	    #plt.figure(figsize=(10,7))
	    
	    #TimeCutoff = 41.0 #kypb
	    #idx_LGM = (rslModelTime <= TimeCutoff)
	    #plt.plot(-rslModelTime[idx_LGM], rslModelAtSite[idx_LGM], '-k', linewidth='1.5')

	    #plt.errorbar(-ageAtSite, rslDataAtSite, yerr=rslerrAtSite, fmt='.r', marker='.', mfc='r', mec='r', ecolor='r')

	    #plt.xlim([-20, 0])
	    #plt.xticks(np.arange(-20,0.3,2.5),['-20','','-15','','-10','','-5','','0'])
	    #plt.margins(0.05,0.05)
	    #plt.title(sitenameAtSite, fontsize=12)
	    #plt.xlabel('Time (kybp)',fontsize=12)
	    #plt.ylabel('RSL (meters)',fontsize=12)

	    #plt.show()
	    #plt.close()

	    # Calculate the RMS between rsl data and model predictions
	    #rslDataError = np.stack((np.abs(rslDataAtSite-rslDataMinAtSite), np.abs(rslDataMaxAtSite-rslDataAtSite)), axis=0).max(axis=0)
        rslDataError = np.abs(rslerrAtSite)
	    #timeDataError = np.stack((np.abs(ageAtSite-ageMinAtSite), np.abs(ageMaxAtSite-ageAtSite)), axis=0).max(axis=0)

	    # Search the index for the nearest model predictions to the rsl data at each time
	    #distanceTime = np.zeros(len(rslModelTime))
	    #idx_rslModel = np.zeros(len(ageAtSite))
	    #time_count = 0

	    #for idx_timeData in ageAtSite:
		#distanceTime = np.abs(rslModelTime-idx_timeData)
		#idx_rslModel[time_count] = distanceTime.argmin()
		#time_count += 1

	    #idx_rslModel = idx_rslModel.astype(int)

	    # Replace the rsl data error with 5.0 meters if the err is less than 5.0
	    # To avoid devided by ~0
	    #rslDataError[rslDataError<1.0] = 1.0

	    # Calculate RMS
	    #err[idx_site] = np.mean(((rslDataAtSite - rslModelAtSite[idx_rslModel]) / rslDataError)**2)

            # Interpolate model predictions at epochs of RSL data point
        f = interpolate.interp1d(rslModelTime, rslModelAtSite)
        rslModel_interp = f(ageAtSite)

            # Compute misfit
        err[idx_site] = np.mean(((rslDataAtSite - rslModel_interp) / rslDataError)**2)
            #pdb.set_trace()
	    # Jacobi


	# Calculate the global average RMS
	errGlobalMean = np.mean(err)
    print('grid {} RMS: {}'.format(idGridFile, errGlobalMean))
	return errGlobalMean

# Generate index for grided files
idGridFileArr = np.arange(1, 865, 1).astype(str)
rms = np.zeros(len(idGridFileArr))
for idx_grid in range(len(idGridFileArr)):
        idGridFile = idGridFileArr[idx_grid].zfill(3)
        
        # Call the function
        rms[idx_grid] = calc_RMS_grid(idGridFile=idGridFile)

#np.savetxt("ICE6G_rms_lambeck_NA_all.txt", rms)
np.savetxt("ICE6G_rms_farfield.txt", rms)
