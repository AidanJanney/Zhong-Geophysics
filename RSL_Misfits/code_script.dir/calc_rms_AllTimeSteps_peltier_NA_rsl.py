#!/usr/bin/evn python
# This code is to read RSL data from Peltier's RSL database
# Read the gridded model predictions
# And plot the curve together to compare
# Then calculate the global averaged RMS

from os import path
import numpy as np
#import matplotlib.pyplot as plt
import pdb
import pandas

def calc_RMS_grid(idGridFile):
	''' 
        calcuale RMS for each grided file
        '''
	# Read RSL gridded model prediction
	#idGridFile = '001'
	#ddirGridFile = '/home/khuan/ICE6G_results_0.1grids/grid'+idGridFile   # ICE6G + 2 layer viscosity
	#ddirGridFile = '/home/khuan/ICE6G_m1_results_0.1grids/grid'+idGridFile+"_m1" # ICE6G_m1
	ddirGridFile = '/home/khuan/ICE6G_vm5a_results/grid'+idGridFile  # ICE6G + vm5a
	#fnGridFile = 'RSL_ICE6G_grids_compr_peltier_rsl.txt'
	#fnGridFile = 'RSL_ICE6G_grids_compr_peltier_digitized.txt'
        fnGridFile = "RSL_ICE6G_grids_compr_peltier_digitized_NA.txt"
        #fnGridFile = "RSL_ICE6G_grids_compr_peltier_digitized_FN.txt"
        #fnGridFile = 'RSL_ICE6G_grids_compr_peltier_NA_ALL.txt'
	rslModel = np.loadtxt(path.join(ddirGridFile, fnGridFile))

	# Read time tag for RSL gridded model prediction
	ddirTime = '/home/khuan/ICE6G_results_0.1grids/grid001'
	fnTime = 'years.reverse'
	rslModelTime = np.loadtxt(path.join(ddirTime, fnTime))

	# Read RSL database from Peltier based on selected UTKEY from file 'UTKEY_pelier_selected.txt'
	# In the selected UTKEY file, the site which is less than 5 data points has been deleted.
	ddirRSL = '/home/khuan/rsl_database/'
	#fnRSL = 'rsl_dbase_peltier_all.xlsx'
        #fnRSL = 'rsl_Peltier_2015paper.xls'
        fnRSL = 'rsl_Peltier_2015paper_NA.xls'
        #fnRSL = 'rsl_Peltier_2015paper_FN.xls'

	#dataPeltier = pandas.read_excel(path.join(ddirRSL, fnRSL), sheet_name='dbase_table')
	dataPeltier = pandas.read_excel(path.join(ddirRSL, fnRSL), sheet_name='sheet 1')
	UTKEY = dataPeltier['UTKEY'].values
	#TYP = dataPeltier['TYP'].values
	RSL = dataPeltier['RSL'].values
	RSLMIN = dataPeltier['RSLMIN'].values
	RSLMAX = dataPeltier['RSLMAX'].values
	AGE = dataPeltier['AGE'].values
	AGEMIN = dataPeltier['AGEMIN'].values
	AGEMAX = dataPeltier['AGEMAX'].values
	LAT = dataPeltier['LAT'].values
	LON = dataPeltier['LON'].values

        data_locations = np.loadtxt("/home/khuan/gridresearch_results_ICE6G_batch/RSL_site_location_Peltier_digitized_NA.txt")
        #data_locations = np.loadtxt("/home/khuan/gridresearch_results_ICE6G_batch/RSL_site_location_Peltier_digitized_FN.txt")
        lats = data_locations[:,1]
        lons = data_locations[:,0]  #Lambeck's longitude range is 0 to 360.

	# Read selected UTKEY file
	ddirUTKEY = '..'
	#fnUTKEY = 'UTKEY_pelier_selected.txt'
        fnUTKEY = 'Sitename_Peltier_digitized_NA.txt'
        #fnUTKEY = 'Sitename_Peltier_digitized_FN.txt'
        #fnUTKEY = 'Sitename_Peltier_NA_all.txt'
	#fnUTKEY = 'UTKEY_pelier_digitized_selected.txt'
	#selectedUTKEY = np.loadtxt(path.join(ddirUTKEY, fnUTKEY))
        f_sitename = open(path.join(ddirUTKEY, fnUTKEY), 'r')
        selectedUTKEY = f_sitename.readlines()

	#err = np.zeros(len(selectedUTKEY))
        #err = np.zeros((len(selectedUTKEY), 3))  # Format: Longitude, Latitude, misfit
        err = np.zeros((len(RSL), 8))  # Format: Longitude, Latitude, time, misfit, prediction, rsl, rslerr, timeerr
        err_idx = 0

	for idx_site in range(len(selectedUTKEY)):
	    #print(idx_site)
            #pdb.set_trace()
	    # Read model predictions
	    rslModelAtSite = rslModel[idx_site, :]
	    rslModelAtSite = rslModelAtSite - rslModelAtSite[-1]

	    # Search rsl data for idx_site
            UTKEYAtSite = selectedUTKEY[idx_site].strip() # delete /n in string
            #UTKEYAtSite = selectedUTKEY[idx_site]
	    idx_UTKEY = (UTKEYAtSite == UTKEY)

            rslDataLon = np.mean(LON[idx_UTKEY])
            rslDataLat = np.mean(LAT[idx_UTKEY])
	    rslDataAtSite = RSL[idx_UTKEY]
	    rslDataMinAtSite = RSLMIN[idx_UTKEY]
	    rslDataMaxAtSite = RSLMAX[idx_UTKEY]

            # ONLY applied Peltier's global rsl database
	    #ageAtSite = AGE[idx_UTKEY]/1000.
	    #ageMinAtSite = AGEMIN[idx_UTKEY]/1000.
	    #ageMaxAtSite = AGEMAX[idx_UTKEY]/1000.

            # Applied on Peltier's digitized database and Lambeck's NA & farfield rsl database
            ageAtSite = AGE[idx_UTKEY]
            ageMinAtSite = AGEMIN[idx_UTKEY]
            ageMaxAtSite = AGEMAX[idx_UTKEY]

	    # Make the plot for comparisons
	    #plt.figure(figsize=(10,7))
	    
	    #TimeCutoff = 20.0 #kypb
	    #idx_LGM = (rslModelTime <= TimeCutoff)
	    #plt.plot(-rslModelTime[idx_LGM], rslModelAtSite[idx_LGM], '-k', linewidth='1.5')

	    #plt.errorbar(-ageAtSite, rslDataAtSite, xerr=np.array([ageMaxAtSite-ageAtSite, ageAtSite-ageMinAtSite]), yerr=np.array([rslDataAtSite-rslDataMinAtSite, rslDataMaxAtSite-rslDataAtSite]), fmt='.r', marker='.', mfc='r', mec='r', ecolor='r')

	    #plt.xlim([-16, 0])
	    #plt.xticks(np.arange(-20,0.3,2.5),['-20','','-15','','-10','','-5','','0'])
	    #plt.margins(0.05,0.05)
	    #plt.title(UTKEYAtSite, fontsize=12)
	    #plt.xlabel('Time (kybp)',fontsize=12)
	    #plt.ylabel('RSL (meters)',fontsize=12)

	    #plt.show()
	    #plt.close()

	    # Calculate the RMS between rsl data and model predictions
	    rslDataError = np.stack((np.abs(rslDataAtSite-rslDataMinAtSite), np.abs(rslDataMaxAtSite-rslDataAtSite)), axis=0).max(axis=0)
	    timeDataError = np.stack((np.abs(ageAtSite-ageMinAtSite), np.abs(ageMaxAtSite-ageAtSite)), axis=0).max(axis=0)

	    # Search the index for the nearest model predictions to the rsl data at each time
	    distanceTime = np.zeros(len(rslModelTime))
	    idx_rslModel = np.zeros(len(ageAtSite))
	    time_count = 0

	    for idx_timeData in ageAtSite:
		distanceTime = np.abs(rslModelTime-idx_timeData)
		idx_rslModel[time_count] = distanceTime.argmin()
		time_count += 1

	    idx_rslModel = idx_rslModel.astype(int)

	    # Replace the rsl data error with 5.0 meters if the err is less than 5.0
	    # To avoid devided by ~0
	    #rslDataError[rslDataError<5.0] = 50.0
            if 0:  
              for i_err in range(len(rslDataError)):
                    if rslDataError[i_err] > 1.0:
                            continue
                    else:
                            if np.abs(rslDataAtSite[i_err]) >= 10.:
                                    rslDataError[i_err] = np.abs(rslDataAtSite[i_err]) * 0.5
                            if np.abs(rslDataAtSite[i_err]) < 10. and np.abs(rslDataAtSite[i_err]) > 1.:
                                    rslDataError[i_err] = np.abs(rslDataAtSite[i_err]) * 1.5
                            if np.abs(rslDataAtSite[i_err]) <= 1.:
                                    rslDataError[i_err] = 10.
	    # Calculate RMS
	    #err[idx_site] = np.mean(((rslDataAtSite - rslModelAtSite[idx_rslModel]) / rslDataError)**2+ ((ageAtSite - rslModelTime[idx_rslModel]) / timeDataError)**2)

	    # Jacobi

            # Compute mean misfit at each site
            #err[idx_site,0] = np.round(rslDataLon,2)
            #err[idx_site,1] = np.round(rslDataLat,2)
            #err[idx_site,2] = np.mean(((rslDataAtSite - rslModelAtSite[idx_rslModel]) / rslDataError)**2+ ((ageAtSite - rslModelTime[idx_rslModel]) / timeDataError)**2)

            # Compute mean misfit at each site
            err[err_idx:err_idx+len(rslDataAtSite),0] = np.round(rslDataLon,2)
            err[err_idx:err_idx+len(rslDataAtSite),1] = np.round(rslDataLat,2)
            err[err_idx:err_idx+len(rslDataAtSite),2] = ageAtSite
            err[err_idx:err_idx+len(rslDataAtSite),3] = ((rslDataAtSite - rslModelAtSite[idx_rslModel]) / rslDataError)**2+ ((ageAtSite - rslModelTime[idx_rslModel]) / timeDataError)**2
            err[err_idx:err_idx+len(rslDataAtSite),4] = rslModelAtSite[idx_rslModel]
            err[err_idx:err_idx+len(rslDataAtSite),5] = rslDataAtSite
            err[err_idx:err_idx+len(rslDataAtSite),6] = rslDataError
            err[err_idx:err_idx+len(rslDataAtSite),7] = timeDataError

            err_idx = err_idx+len(rslDataAtSite)

        # Return RMS dictribution at RSL sites for the input model
        return err



	# Calculate the global average RMS
	#errGlobalMean = np.mean(err)
	#print('grid {} RMS: {}'.format(idGridFile, errGlobalMean))
        #return errGlobalMean

# Generate index for grided files
#idGridFileArr = np.arange(1, 865, 1).astype(str)
#rms = np.zeros(len(idGridFileArr))
#for idx_grid in range(len(idGridFileArr)):
#        idGridFile = idGridFileArr[idx_grid].zfill(3)
        
#        # Call the function
#        rms[idx_grid] = calc_RMS_grid(idGridFile=idGridFile)

#np.savetxt("ANU_rms_Peltier_rsl.txt", rms)
#np.savetxt("ANU_rms_Peltier_rsl_digitized_NA_timeErr.txt", rms)
#np.savetxt("ANU_rms_Peltier_rsl_NA_ALL.txt", rms)

# Calculate the rms distribution at each rsl sites for the best fit model grids563 for ANU + Lameck's NA rsl dataset
#rms = calc_RMS_grid(idGridFile="565")
#rms = calc_RMS_grid(idGridFile="661")
#rms = calc_RMS_grid(idGridFile="594") # ICE6G_m1
rms = calc_RMS_grid(idGridFile="001") # ICE6G_vm5a

#np.savetxt("ICE6G_rms_presentday_Peltier_digitized_NA_594.txt", rms, fmt="%0.2f  %0.2f  %6.3f")
#np.savetxt("ICE6G_m1_rms_alltimesteps_Peltier_digitized_FN_594.txt", rms, fmt="%0.2f  %0.2f  %6.3f  %6.3f     %6.3f  %6.3f   %6.3f  %6.3f")
np.savetxt("ICE6G_vm5a_rms_alltimesteps_Peltier_digitized_NA.txt", rms, fmt="%0.2f  %0.2f  %6.3f  %6.3f     %6.3f  %6.3f   %6.3f  %6.3f")
#np.savetxt("ICE6G_rms_alltimesteps_Peltier_digitized_NA_594.txt", rms, fmt="%0.2f  %0.2f  %6.3f  %6.3f     %6.3f  %6.3f   %6.3f  %6.3f")
#np.savetxt("ICE6G_rms_alltimesteps_Peltier_digitized_FN_594.txt", rms, fmt="%0.2f  %0.2f  %6.3f  %6.3f     %6.3f  %6.3f   %6.3f  %6.3f")
