#!/usr/bin/evn python
# This code is to read RSL data from Lambeck's RSL database
# Read the gridded model predictions
# And plot the curve together to compare
# Then calculate the global averaged RMS

from os import path
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pdb
import pandas
from matplotlib.ticker import MultipleLocator

def calc_RMS_grid(idGridFile):
 ''' 
        calcuale RMS for each grided file
        '''
 # Read RSL gridded model prediction
 #idGridFile = '001'
 #ddirGridFile = f'/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/A{a:.2f}L{l:.2f}/'+idGridFile
 ddirGridFile = f'/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/U{u:.2f}L{l:.2f}/'+idGridFile
 #fnGridFile = f"/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/A{a:.2f}L{l:.2f}/RSL_ICE6G_grids_compr_lambeck_farfield.txt"
 fnGridFile = f"/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/U{u:.2f}L{l:.2f}/RSL_ICE6G_grids_compr_lambeck_farfield.txt"
 #fnGridFile = f"/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/U{u:.2f}L{l:.2f}/RSL_ICE6G_grids_compr_lambeck_NA.txt"
 #fnGridFile = 'RSL_ICE6G_grids_compr_lambeck_NA_ALL.txt'
 rslModel = np.loadtxt(fnGridFile)

 # Read time tag for RSL gridded model prediction
 #ddirTime = '/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/A18.50L21.00'
 ddirTime = '/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/U19.00L21.00'
 fnTime = 'years.reverse'
 rslModelTime = np.loadtxt(path.join(ddirTime, fnTime))

 # Read RSL database from Lambeck based on selected sitename from file 'Sitename_Lambeck_FarField_selected.txt'
 # In the selected sitename file, the site which is less than 5 data points has been deleted.
 ddirRSL = '/home/sebi8707/Analytic_GIA_Code/version3/Prj_Sedi/RSL_OBS/RSL_dataset_Lambeck_farfield'
 #ddirRSL = '/home/sebi8707/Analytic_GIA_Code/version3/Prj_Sedi/RSL_OBS/RSL_dataset_Lambeck_NA'
 
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
     #plt.savefig(f'/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/RSLmisfit_{idx_site}.png')

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

#Save data to file
rms_values = []
with open('output_data.txt', 'w') as file:
    for l in [round(x, 2) for x in np.arange(21.00, 23.01, 0.25)]:
        for u in [round(x, 2) for x in np.arange(19.00, 21.01, 0.25)]:
        #for a in [round(x, 2) for x in np.arange(18.50,21.01, 0.25)]:
            #idGridFile = f'A{a:.2f}L{l:.2f}'
            idGridFile = f'U{u:.2f}L{l:.2f}'
            rms_value = calc_RMS_grid(idGridFile=idGridFile)
            rms_values.append(rms_value)
            #file.write(f'{a}\t{l}\t{rms_value}\n')
            file.write(f'{u}\t{l}\t{rms_value}\n')

# Convert the list of RMS values to a NumPy array
rms = np.array(rms_values)

#plot
data = np.loadtxt('output_data.txt')
#asth_visc = data[:, 0]
um_visc = data[:, 0]
lm_visc = data[:, 1]
rms = data[:, 2]

#reshape rms
lm_values = np.arange(21.00, 23.01, 0.25)
#a_values = np.arange(18.50,21.01, 0.25)
um_values = np.arange(19.00,21.01, 0.25)
#rms_2D = np.reshape(rms, (len(lm_values),len(a_values)))
rms_2D = np.reshape(rms, (len(lm_values),len(um_values)))
transposed_rms = list(map(list, zip(*rms_2D)))

#fig, ax = plt.subplots(figsize=(12, 8))
fig, ax = plt.subplots(figsize=(14, 10))

#scatter = ax.scatter(lm_visc, asth_visc, c=rms, cmap='viridis', marker='o', s=50)
#cbar = plt.colorbar(scatter)
#cbar.set_label('RMS Error')

#misfit_plot = plt.imshow(rms_2D, vmin=0, vmax=60, cmap='seismic')
misfit_plot = plt.imshow(transposed_rms, vmin=10, vmax=70, cmap='jet')
cbar = plt.colorbar(misfit_plot)
cbar.set_label('RMS Error')
cbar.ax.yaxis.set_major_locator(MultipleLocator(5))

#Set axes ticks
#ax.set_yticks(np.arange(len(np.unique(lm_visc))))
#ax.set_xticks(np.arange(len(np.unique(asth_visc))))
ax.set_xticks(np.arange(len(np.unique(lm_visc))))
ax.set_yticks(np.arange(len(np.unique(um_visc))))
#ax.set_yticklabels(np.unique(lm_visc))
#ax.set_xticklabels(np.unique(asth_visc))
ax.set_xticklabels(np.unique(lm_visc))
ax.set_yticklabels(np.unique(um_visc))
plt.gca().invert_yaxis()

# Set labels and title
#ax.set_ylabel('Lower Mantle Viscosity')
#ax.set_xlabel('Asthenosphere Viscosity')
ax.set_xlabel('Lower Mantle Viscosity')
ax.set_ylabel('Upper Mantle Viscosity')
#plt.title()

# Show legend
#ax.legend()

# Show the plot
#plt.show()
#plt.savefig(f'/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/RSLmisfit.png')
plt.savefig(f'/home/sebi8707/Analytic_GIA_Code/version3/Prj_Te_WeakAstheno/copy_Multi_Analytic_cases/cases/UMLM_RSLmisfit_LambeckFarfield.png')

#Print rms
#print(rms)

#np.savetxt("ICE6G_rms_lambeck_NA_all.txt", rms)
np.savetxt("ICE6G_rms_farfield.txt", rms)
