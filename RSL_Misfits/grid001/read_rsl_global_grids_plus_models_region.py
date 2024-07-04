#!/usr/bin/evn python
#this code is to read Lambeck's rsl and Peltier's excel table for 1 by 1 global gird 
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas
import matplotlib
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

ddir_of = './' #'./fig_rsl_global_grids_all_region_shorter'

#generate the global grids
lats = np.arange(90.,-91.,-1.0) #np.arange(-90.0,91.0,1.0)
lons = np.arange(0.,361.,1.0)

#read the peltier's table
data1 = pandas.read_excel('../../rsl_dbase_peltier_all.xlsx',sheet_name = 'dbase_table')
UTKEY1 = data1['UTKEY'].values
TYP1 = data1['TYP'].values
RSL1 = data1['RSL'].values
RSLMIN1 = data1['RSLMIN'].values
RSLMAX1 = data1['RSLMAX'].values
AGE1 = data1['AGE'].values
AGEMIN1 = data1['AGEMIN'].values
AGEMAX1 = data1['AGEMAX'].values
LAT1 = data1['LAT'].values
LON1 = data1['LON'].values
idx = LON1 < 0.0
LON1[idx] = LON1[idx] + 360.0   #longitude range is -180 to 180 in Peltier's table

#read the peltier's table from digitized rsl in his2015 paper
data3 = pandas.read_excel('../../RSL_measurements/rsl_Peltier_2015paper.xls',sheet_name = 'sheet 1')
UTKEY3 = data3['UTKEY'].values
RSL3 = data3['RSL'].values
RSLMIN3 = data3['RSLMIN'].values
RSLMAX3 = data3['RSLMAX'].values
AGE3 = data3['AGE'].values
AGEMIN3 = data3['AGEMIN'].values
AGEMAX3 = data3['AGEMAX'].values
LAT3 = data3['LAT'].values
LON3 = data3['LON'].values

#read the lambeck's farfield table(from his 2014 paper)
data = pandas.read_excel('../../Lambeck_farfield_RSL/Lambeck_etal_2014_PNAS_supp_edit.xlsx')
SITENAME = data['Site'].values
RSL = data['RSL'].values
RSLERR = data['Sig.RSL'].values
AGE = data['Age'].values
LAT = data['Latitude'].values
LON = data['Longitude'].values
SAMPLE = data['Sample'].values

#read the lambeck's NA table (from his 2017 paper)
data2 = pandas.read_excel(path.join('../../ANU_ice_model/QSR_SI/rsl_observations','rsl data set (after Dyke 2004).xlsx'))
SITENAME2 = data2['site name'].values
RSL2 = data2['observed rsl'].values
RSLERR2 = data2['2xsigma rsl'].values
AGE2 = data2['time ka BP'].values
LAT2 = data2['latitude N'].values
LON2 = data2['longitude E'].values

#idx = LON >= 180.0
#LON[idx] = LON[idx] - 360.0


#plot the map of the rsl sites from different groups
plt.figure()
plt.title('''RSL sites: Peltier_all(red); Lambeck_NA(green); 
Lambeck_Far(blue); Peltier_digitized(yellow)''')
m = Basemap()
m.drawcoastlines(zorder=1)
#pdb.set_trace()
#Peltier's rsl sites
plt.scatter(LON1[LON1>=180.]-360.,LAT1[LON1>=180.],s=10,marker='.',color='red')
plt.scatter(LON1[LON1<180.],LAT1[LON1<180.],s=10,marker='.',color='red')

#Lameck's far field rsl sites
plt.scatter(LON[LON>=180.]-360.,LAT[LON>=180.],s=10,marker='.',color='blue')
plt.scatter(LON[LON<180.],LAT[LON<180.],s=10,marker='.',color='blue')

#Lambeck's NA rsl sites
plt.scatter(LON2[LON2>=180.]-360.,LAT2[LON2>=180.],s=10,marker='.',color='green')
plt.scatter(LON2[LON2<180.],LAT2[LON2<180.],s=10,marker='.',color='green')

#Peltier's digitized rsl sites from his 2015 paper
plt.scatter(LON3[LON3>=180.]-360.,LAT3[LON3>=180.],s=10,marker='x',color='yellow')
plt.scatter(LON3[LON3<180.],LAT3[LON3<180.],s=10,marker='x',color='yellow')

#plt.savefig(path.join(ddir_of,'map_all_rsl_site.png'),bbox_inches='tight',pad_inches = 0)
#plt.show()
plt.close()

#limits for lat and lon centered at each grid
lim_lon = 0.5 #1.0 #0.5
lim_lat = 0.5 #1.0 #0.5

# txt for record the rsl location
#region = 'Farfield' #NA 1; GreenLand 2; FN 3; Antarctica 4; Farfield 0
#region_id = 0 
#fn_txt = open(path.join(ddir_of,'rsl_location_record_predictions_'+region+'.txt'),'w')

#open region mask
#region_mask = np.loadtxt('../ice6g_data/region_tag.txt')

#read model prediction case60t
rsl_case60t = np.loadtxt('./RSL_ICE6G_vm5a_compr.txt')
time_case60t = np.loadtxt('./years.reverse')
# rsl relative to present day
repetitions = len(rsl_case60t[0])
rsl_case60t = rsl_case60t - np.transpose([rsl_case60t[:,-1]]*repetitions)


#initial error array
err = np.zeros(len(lats)*len(lons))
err_Jacobi = np.zeros(len(lats)*len(lons))
Jacobi = np.zeros(len(lats)*len(lons))

# read ice mask at LGM to seperate the sites in the near and far field
ice_LGM = np.loadtxt('../../ice6g_data/ICE1x1_26ka/ice6g180x360.1')

#for each grid search the data around the grid 1 by 1 degree
count = -1
count_grid = -1

for lat in lats:
	for lon in lons:
#test for tahiti
#for lat in [-17.,-16.,]:
#	for lon in [210.,209.]:
		count_grid += 1
                count += 1

		#if region_mask[count,2] != region_id:
		#	continue;

		#group peltier's table
		idx_lon1 = (LON1 >= lon - lim_lon) & (LON1 <= lon + lim_lon)
		idx_lat1 = (LAT1 >= lat - lim_lat) & (LAT1 <= lat + lim_lat)
		idx1 = idx_lon1 & idx_lat1

		site_rsl = RSL1[idx1]
		site_rsl_min = RSLMIN1[idx1]
		site_rsl_max = RSLMAX1[idx1]
		site_age = AGE1[idx1]/1000.
		site_age_min = AGEMIN1[idx1]/1000.
		site_age_max = AGEMAX1[idx1]/1000.
		site_UTKEY = UTKEY1[idx1]
		#print('Peltier index')
		#print(site_UTKEY)

               #group peltier's digitized rsl table
                idx_lon3 = (LON3 >= lon - lim_lon) & (LON3 <= lon + lim_lon)
                idx_lat3 = (LAT3 >= lat - lim_lat) & (LAT3 <= lat + lim_lat)
                idx3 = idx_lon3 & idx_lat3

                site_rsl3 = RSL3[idx3]
                site_rsl_min3 = RSLMIN3[idx3]
                site_rsl_max3 = RSLMAX3[idx3]
                site_age3 = AGE3[idx3]
                site_age_min3 = AGEMIN3[idx3]
                site_age_max3 = AGEMAX3[idx3]
                site_UTKEY3 = UTKEY3[idx3]
                #print('Peltier digitied rsl index')
                #print(site_UTKEY3)

		#group lambeck's farfield table
                idx_lon = (LON >= lon - lim_lon) & (LON <= lon + lim_lon)
                idx_lat = (LAT >= lat - lim_lat) & (LAT <= lat + lim_lat)
                idx = idx_lon & idx_lat

		sitename_lambeck = SITENAME[idx]
		rsl_lambeck = RSL[idx]
		rslerr_lambeck = RSLERR[idx]
		age_lambeck = AGE[idx]
		lat_lambeck = LAT[idx]
		lon_lambeck = LON[idx]
		sample_lambeck = SAMPLE[idx]
		#print('''Lambeck's far field''')
		#print(sitename_lambeck)

                #group lambeck's NA table
                idx_lon2 = (LON2 >= lon - lim_lon) & (LON2 <= lon + lim_lon)
                idx_lat2 = (LAT2 >= lat - lim_lat) & (LAT2 <= lat + lim_lat)
                idx2 = idx_lon2 & idx_lat2

                sitename_lambeck2 = SITENAME2[idx2]
                rsl_lambeck2 = RSL2[idx2]
                rslerr_lambeck2 = RSLERR2[idx2]
                age_lambeck2 = AGE2[idx2]
                lat_lambeck2 = LAT2[idx2]
                lon_lambeck2 = LON2[idx2]
                #print('''Lambeck's NA rsl''')
                #print(sitename_lambeck2)

		#LGM ice
               	idx_lon_ice = (ice_LGM[:,0] >= lon - lim_lon) & (ice_LGM[:,0] <= lon + lim_lon)
                idx_lat_ice = (ice_LGM[:,1] >= lat - lim_lat) & (ice_LGM[:,1] <= lat + lim_lat)
                idx_ice = idx_lon_ice & idx_lat_ice
		ih_select = ice_LGM[idx_ice,2]
		ih_ave = np.mean(ih_select)

		#select near or far field sites by ice height 
		#if ih_ave >= 1.0:   # ih_ave < 1.0 meter:  for near field
		#	continue

		#make plots
		if len(site_rsl)==0 and len(rsl_lambeck)==0 and len(rsl_lambeck2)==0:
			continue
		#pdb.set_trace()
		if len(site_rsl)==0 and len(rsl_lambeck)!=0:
			if np.max(age_lambeck) < 6.0:
				continue

			if len(age_lambeck) <= 8:
				continue

                if len(site_rsl)!=0 and len(rsl_lambeck)==0:
	                if np.max(site_age) < 6.0:
        	                continue

                	if len(site_age) <= 8:
                        	continue

                if len(site_rsl)!=0 and len(rsl_lambeck)!=0:
                	if np.max(site_age) < 6.0 and np.max(age_lambeck) < 6.0:
                        	continue

                	if len(site_age) <= 8 and len(age_lambeck) <= 8:
                        	continue

		print(count)
		#plot the map of the data coordinate for i
		#pdb.set_trace()
		fig = plt.figure(figsize=(10,5)) #,constrained_layout=True)
		#fig.set_figheight(6)
		#fig.set_figwidth(8)
		#ax1 = plt.subplot2grid(shape=(3,5),loc=(0,3),colspan=2,rowspan=2)
		ax1 = fig.add_axes([0.65,0.38,0.3,0.52])
		#gs = fig.add_gridspec(3,4)
		#ax1 = fig.add_subplot(gs[0,-1])
		#axs,(ax1,ax2) = fig.subplots(1,2,gridspec_kw={'width_ratios':[2, 1],'height_ratios':[2,1]})
		#ax1.set_title('site location')
		m = Basemap()

		#plot the ice mask at LGM
		ih = ice_LGM[:,2]
		idx = (ih > 1.0)  #1.0 meter
		ih[idx] = 1.0 # marked the region of ice 
		ih=ih.reshape([180,360])
		shp = np.shape(ih)
		ih1 = np.zeros((shp[0],shp[1]))

                lons_ice = np.arange(0.5,360.,1.)
                lats_ice = np.arange(89.5,-90,-1.0)	

		#lons = lons-180.0
		for i_lon in range(len(lons_ice)):
		        if i_lon <= 180:
                		ih1[:,i_lon] = ih[:,i_lon+179]
		for i_lon in range(len(lons_ice)):
      			if i_lon > 180:
                		ih1[:,i_lon] = ih[:,i_lon-180]

		
		cmap1 = LinearSegmentedColormap.from_list("my_colormap", ((1, 1, 1), (0.6, 0.6, 0.6)), N=6, gamma=1.0)
		ax1.pcolor(lons_ice-180.0,lats_ice,ih1,vmin=0,vmax=1.,cmap=cmap1,zorder=1) 
		#ax1.pcolor(lons-180.0,lats,ih1,cmap=cmap1,zorder=1)
                m.drawcoastlines(color='black',zorder=2)	

		#plot the site location marked by x
		#xi,yi = m(lon,lat)
		if lon >= 180.0:
			lon = lon - 360.0
		ax1.scatter(lon,lat,s=30,marker='x',color='red',alpha=1.0,zorder=3)
		#m.drawcoastlines()
		#plt.show()
		#pdb.set_trace()
		#if (np.max(site_age) >= time_cutoff) or (np.max(AGE[idx_site]) >= time_cutoff):
		#	plt.savefig(path.join(ddir_of,'FAR'+str(count)+'_site_position.png'),bbox_inches='tight',pad_inches = 0)
		#plt.close()


		#plt.figure()
		#model predictions
		#ax2 = fig.add_subplot(gs[:,:-1])
		#ax2 = plt.subplot2grid(shape=(3,5),loc=(0,0),colspan=3,rowspan=3)
	        ax2 = fig.add_axes([0.1,0.11,0.53,0.79])

		#ice6g+vm5a
		cutoff_time = 16.0 #21.0 #26.0
        	idx_LGM = time_case60t <= cutoff_time
        	ax2.plot(-1*time_case60t[idx_LGM],rsl_case60t[count_grid,idx_LGM],'-k',linewidth=1.5)

		#peltier rsl
		ax2.errorbar(-site_age,site_rsl,xerr=np.array([site_age_max-site_age,site_age-site_age_min]),yerr=np.array([site_rsl-site_rsl_min,site_rsl_max-site_rsl]),fmt='.r', marker='.', mfc='r',mec='r', ecolor='r')
		
		#Lambeck rsl
		ax2.errorbar(-age_lambeck,rsl_lambeck,xerr=0.,yerr=rslerr_lambeck,fmt='.b', marker='.', mfc='b',mec='b', ecolor='b')
		
                #Lambeck rsl
                ax2.errorbar(-age_lambeck2,rsl_lambeck2,xerr=0.,yerr=rslerr_lambeck2,fmt='.g', marker='.', mfc='g',mec='g', ecolor='g')

                #peltier digitized rsl
                #ax2.errorbar(-site_age3,site_rsl3,xerr=np.array([site_age_max3-site_age3,site_age3-site_age_min3]),yerr=np.array([site_rsl3-site_rsl_min3,site_rsl_max3-site_rsl3]),fmt='.y', marker='.', mfc='y',mec='y', ecolor='y')

		ax2.set_xlabel('Time (kybp)')
		ax2.set_ylabel('RSL (meters)')
		ax2.set_title(str(count))
		ax2.set_xlim([-15.,0.])
		#plt.ylim([-150.,280.])  #([-310.,280.])
		#plt.savefig(path.join(ddir_of,'RSL_peltier_lambeck_'+str(count).zfill(5)+'_predictions_map_'+region+'.png'))  #,bbox_inches='tight',pad_inches = 0)
		#plt.show()
		plt.close()

		#compute the rms of data-model misfit
		#combine the RSL from Lambeck and Peltier
		#pdb.set_trace()
		rsl_all = np.concatenate((site_rsl,rsl_lambeck,rsl_lambeck2),axis=0)
		rsl_time_all = np.concatenate((site_age,age_lambeck,age_lambeck2),axis=0)
		rsl_err_all = np.concatenate((np.stack((np.abs(site_rsl-site_rsl_min),np.abs(site_rsl_max-site_rsl)),axis=0).max(axis=0),rslerr_lambeck,rslerr_lambeck2),axis=0)

	        #time for rsl data      
        	time_R = -1*rsl_time_all

		#time for case60t
		time_A = -1*time_case60t
		diff_t = np.zeros(len(time_A))
		idx_R = np.zeros(len(time_R))
		#pdb.set_trace()
		for i_time in range(len(time_R)):
			diff_t = np.abs(time_A - time_R[i_time])
			idx_R[i_time] = diff_t.argmin()

		idx_R_60t = idx_R.astype(int)

		idx_rsl_err = rsl_err_all < 5.0 
		rsl_err_all[idx_rsl_err] = 5.0

		#pdb.set_trace()
      		err[count] = np.mean(((rsl_all - rsl_case60t[count_grid,idx_R_60t])/rsl_err_all)**2)
		Jacobi[count] = np.sin((90.-lat)*np.pi/180.)*1.*np.pi/180.*1.*np.pi/180.
				
	#	plt.text(0.01,0.01,'ANU+NA_visc: '+str(round(err13[count],2))+'; ANU+FAR_visc: '+str(round(err14[count],2))+'; ICE6g+vm5a: '+str(round(err60tc[count],2)),horizontalalignment='left',verticalalignment='bottom', transform=ax.transAxes)
                #str_title = 'ANU+NA_visc: '+str(round(err13[count],2))+'; ANU+FAR_visc: '+str(round(err14[count],2))+'; ICE6g+vm5a: '+str(round(err60tc[count],2))
                str_title = 'Lat: '+str(lat)+' Lon: '+str(lon)
		plt.title('Grid index: '+str(count)+' ( '+str_title+' )')
	#	plt.savefig(path.join(ddir_of,'RSLdata_model_comparison'+str(count).zfill(5)+'_'+region+'_rms.png'))  #,bbox_inches='tight',pad_inches = 0)
                #plt.show()

                plt.close()
		

		#record the rsl in txt
		#pdb.set_trace()
		#fn_txt.write('%s\n'%('-------------'))
		#fn_txt.write('%i\n'%(count))
		#fn_txt.write('%5.2f  %5.2f\n'%(lon,lat))		
		#fn_txt.write('%s\n'%(site_UTKEY))
		#fn_txt.write('%s\n'%(site_UTKEY3))
		#fn_txt.write('%s\n'%(sitename_lambeck))
		#fn_txt.write('%s\n'%(sitename_lambeck2))
		#count += 1

#plot the rsl on the map for each models
#pdb.set_trace()
#err = [err13,err14,err60tc]
err_mean = np.sum(Jacobi*err)/np.sum(Jacobi)
print err_mean
#np.savetxt(path.join(ddir_of,'rms_'+region+'.txt'),err)

