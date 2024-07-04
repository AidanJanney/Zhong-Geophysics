#!/usr/bin/env pyithon
#This code is to plot interpolated  RSL for 1 Newtonian and 3 non-Newtonian cases with different transition stress in casevm5a
#this code is to plot the Fig. 7 and Fig. 9 in the paper
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas
import matplotlib
from scipy import interpolate

matplotlib.rc('xtick',labelsize=12)
matplotlib.rc('ytick',labelsize=12)

#sitename = 'Churchill' #'Sotra_Norway'  #'Barbados' #'Boston' #'Richmond_Gulf' #'Churchill' #'Richmond_Gulf'
#sitename = 'NewYork' #'NewYork' #'Murray_Bay' #'Pt_Caen'

#N_site = 18 #12 #18  #North America 
#Tag_site = 'N' #'F' #'N'
#arr_sitename = ['Tahiti','Kiritimati','Christchurch','Port Pirie','Redcliff Belperio','Port Gawler','Fisherman Bay','Wood Point','Abrolhos','Grub Reef','Pioneer Bay','Magnetic Island','Yule Point','Karumba','Huon Peninsula','Semakau','Geylang','Phuket','Ca_Na','Maldives Maalhosmadulu','Maldives Rasdhoo','Reunion','Seychelles','Mayotte','Natal','Barbados']

#arr_sitename = ['W Ungava Bay','Deception Bay','Kugluktuk','Riviere du Loup','Goulburn Lake','Boston','New York','Frosta Norway','Ski Moraine Norway','Sandsjobacka Sweden','Kragero Norway','Sotra Norway','Borgan Norway','Vestfonna','Barbados','Huon Peninsula','Mayotte']#'Vestfonna','Spitsbergen']

#loc_arr = np.array([[290.0,59.4],[285.6,62.1],[244.9,67.9],[290.5,47.8],[252.0,67.0],[289.0,42.0],[286.1,41.2],[10.9,63.6],[10.8,59.8],[12.0,57.4],[9.4,58.9],[5.1,60.3],[11.0,65.0],[19.0,80.0],[300.45,13.04],[147.62,-6.11],[45.27,-12.80]])

#arr_sitename = ['W Ungava Bay','Deception Bay','Kugluktuk','Riviere du Loup','Goulburn Lake','S Massachusetts','New York','Ft George','Frosta Norway','Ski Moraine Norway','Sandsjobacka Sweden','Kragero Norway','Sotra Norway','Borgan Norway','Vestfonna','Varangerfjord Norway','Barbados','Huon Peninsula','Mayotte']

#loc_arr = np.array([[290.0,59.4],[285.6,62.1],[244.9,67.9],[290.5,47.8],[252.0,67.0],[289.5,41.6],[286.1,41.2],[282.1,53.6],[10.9,63.6],[10.8,59.8],[12.0,57.4],[9.4,58.9],[5.1,60.3],[11.0,65.0],[19.0,80.0],[29.4,70.1],[300.45,13.04],[147.62,-6.11],[45.27,-12.80]])

#arr_sitename = ['N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','Tahiti','Kiritimati','Christchurch','Port Pirie','Redcliff Belperio','Port Gawler','Fisherman Bay','Wood Point','Abrolhos','Grub Reef','Pioneer Bay','Magnetic Island','Yule Point','Karumba','Huon Peninsula','Semakau','Geylang','Phuket','Ca_Na','Maldives Maalhosmadulu','Maldives Rasdhoo','Reunion','Seychelles','Mayotte','Natal','Barbados']
#loc_arr = np.array([[304.2,51.5],[295.7,59.8],[290.0,59.4],[285.6,62.1],[282.1,53.6],[265.6,58.7],[244.9,67.9],[232.0,52.0],[236.8,49.6],[263.2,28.2],[269.8,30.0],[279.8,25.4],[279.6,32.7],[286.1,41.2],[289.5,41.6],[289.9,43.6],[295.7,45.8],[290.5,47.8],[29.4,70.1],[19.0,69.8],[10.9,63.6],[5.1,60.3],[9.4,58.9],[10.8,59.8],[12.0,57.4],[17.7,59.5],[19.9,64.0],[23.3,60.1],[24.9,59.4],[210.45,-17.45],[202.52,1.98],[172.67,-43.5],[138.01,-33.16],[137.87,-32.69],[138.46,-34.64],[137.88,-33.43],[137.86,-33.33],[113.83,-28.68],[147.43,-18.63],[146.48,-18.60],[146.87,-19.15],[145.52,-16.57],[140.83,-17.42],[147.62,-6.11],[103.77,1.21],[103.87,1.31],[98.41,7.75],[108.84,11.33],[73.03,5.27],[72.98,4.30],[55.23,-21.08],[55.52,-4.68],[45.27,-12.80],[32.40,-28.47],[300.45,13.04]])

arr_sitename = ['N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11'] #,'Tahiti','Kiritimati','Christchurch','Port Pirie','Redcliff Belperio','Port Gawler','Fisherman Bay','Wood Point','Abrolhos','Grub Reef','Pioneer Bay','Magnetic Island','Yule Point','Karumba','Huon Peninsula','Semakau','Geylang','Phuket','Ca_Na','Maldives Maalhosmadulu','Maldives Rasdhoo','Reunion','Seychelles','Mayotte','Natal','Barbados']
loc_arr = np.array([[304.2,51.5],[295.7,59.8],[290.0,59.4],[285.6,62.1],[282.1,53.6],[265.6,58.7],[244.9,67.9],[232.0,52.0],[236.8,49.6],[263.2,28.2],[269.8,30.0],[279.8,25.4],[279.6,32.7],[286.1,41.2],[289.5,41.6],[289.9,43.6],[295.7,45.8],[290.5,47.8],[29.4,70.1],[19.0,69.8],[10.9,63.6],[5.1,60.3],[9.4,58.9],[10.8,59.8],[12.0,57.4],[17.7,59.5],[19.9,64.0],[23.3,60.1],[24.9,59.4]]) #,[210.45,-17.45],[202.52,1.98],[172.67,-43.5],[138.01,-33.16],[137.87,-32.69],[138.46,-34.64],[137.88,-33.43],[137.86,-33.33],[113.83,-28.68],[147.43,-18.63],[146.48,-18.60],[146.87,-19.15],[145.52,-16.57],[140.83,-17.42],[147.62,-6.11],[103.77,1.21],[103.87,1.31],[98.41,7.75],[108.84,11.33],[73.03,5.27],[72.98,4.30],[55.23,-21.08],[55.52,-4.68],[45.27,-12.80],[32.40,-28.47],[300.45,13.04]])


err = np.zeros(len(arr_sitename))
err_Jacobi = np.zeros(len(arr_sitename))
Jacobi = np.zeros(len(arr_sitename))

ddir_of = './' #'./fig_proposal2022b'  #'./fig_case70' #'./fig_farfield.dir'

#ddir = '../Lambeck_farfield_RSL'
#fn = 'RSL_farfield.xlsx'
#xls = pandas.ExcelFile(path.join(ddir,fn))
#sheets = xls.sheet_names

#set upper and lower limits for the RSL data errorbar
#North America
#upperlims=[17,0,8,10,4,4,7,27,1,11,27,0,3,19,13,6,23,14]
#lowerlims=[21,3,13,26,11,13,12,1,6,14,25,21,0,5,12,17,11,26]

#Fennoscandia
#upperlims=[2,0,0,31,0,0,0,9,0,13,6,0]
#lowerlims=[1,0,0,7,2,7,0,5,0,1,6,0]

#readd the peltier's table
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

#limits for lat and lon centered at each grid
lim_lon = 0.5 #0.25 #1.0 #0.5
lim_lat = 0.5 #0.25 #1.0 #0.5

for i_site in range(len(arr_sitename)):
#for sitename in ['Barbados','Barbados']:
	sitename = arr_sitename[i_site]
	print(sitename)

	lon = loc_arr[i_site,0]
	lat = loc_arr[i_site,1]

       #read compressible ice6g+vm5a
        ddir_compr= './' #'../ice6g_incompressible_ANU.dir/results_50p12/ts_sites_output.dirr' #'../check_code/case_compr_vm5a'
        caseid= 'casevm5a_compr' #'case_compr_vm5a'
        fn_rsl0 = sitename+'_rsl_'+caseid+'.txt'
        fn_g0 = sitename+'_geoid_'+caseid+'.txt'
        fn_u0 = sitename+'_uplift_'+caseid+'.txt'

        topo_g01=np.loadtxt(path.join(ddir_compr,fn_g0))
        topo_u01=np.loadtxt(path.join(ddir_compr,fn_u0))
        topo_rsl01=np.loadtxt(path.join(ddir_compr,fn_rsl0))
        topo_rsl01[:,1] = topo_rsl01[:,1] - topo_rsl01[-1,1]

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



	#read RSL far field data from Lambeck 2014 paper
        #data = pandas.read_excel(path.join(ddir,fn),sheet_name = sitename)
        #age = data['Age'].values
        #rsl = data['RSL'].values
        #sigma = data['Sig.RSL'].values

	#read peltier 2015 Barbados rsl data
	#rsl_peltier = np.loadtxt('/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/RSL_measurements/Barbados.txt')
        #y_error2 = np.copy(rsl_peltier[:,3])
        #y_error1 = np.zeros(len(rsl_peltier[:,3]))
        #y_error = [y_error1, y_error2]


	plt.figure()
	#pdb.set_trace()

        #-----ICE6G+vm5a  incompressible models with vm5a
        cutoff_time = 20.0
        idx_LGM = np.where(-topo_rsl01[:,0]==cutoff_time)
        plt.plot(topo_rsl01[idx_LGM[0][0]:,0],topo_rsl01[idx_LGM[0][0]:,1],'--k',linewidth=1.5)


       #peltier rsl
	plt.errorbar(-site_age,site_rsl,xerr=np.array([site_age_max-site_age,site_age-site_age_min]),yerr=np.array([site_rsl-site_rsl_min,site_rsl_max-site_rsl]),fmt='.g', marker='.', mfc='g',mec='g', ecolor='g')

	#peltier digitized rsl
	#plt.errorbar(-site_age3,site_rsl3,xerr=np.array([site_age_max3-site_age3,site_age3-site_age_min3]),yerr=np.array([site_rsl3-site_rsl_min3,site_rsl_max3-site_rsl3]),fmt='.y', marker='.', mfc='y',mec='y', ecolor='y')
	
	#Lambeck rsl
	plt.errorbar(-age_lambeck,rsl_lambeck,xerr=0.,yerr=rslerr_lambeck,fmt='.b', marker='.', mfc='b',mec='b', ecolor='b')

	#Lambeck rsl
	plt.errorbar(-age_lambeck2,rsl_lambeck2,xerr=0.,yerr=rslerr_lambeck2,fmt='.r', marker='.', mfc='r',mec='r', ecolor='r')


	#plt.legend(['ICE6g+Non_Nwt','ICE6g+vm5a','ANU+Lambeck_NA','ANU+Lambeck_FF','ICE6g+Lambeck_FF'],frameon=False, prop={'size':12})
	#plt.legend(['ICE6g+Non_Nwt','ICE6g+vm5a','ANU+Lambeck_NA','ANU+Lambeck_FF'],frameon=False, prop={'size':12})
        #plt.legend(['ICE6G+vm5a','ICE6G+Lambeck1','W12+vm5a','ICE7G+vm7','ICE6G+Lambeck2'])
	#plt.legend(['Case 0','Case 1','Case 2','Case 3','Case 4'],frameon=False, prop={'size':20})
	#if i_site == 0:
	#	plt.legend(['case 0','case 1','case 2','case 3'],frameon=False, prop={'size':20})	

        #if Tag_site == 'N':
        #        if i_site == 5 or i_site == 13:
	#        	plt.legend([line_case4],['case 4'],frameon=False, prop={'size':20})
	
	plt.xlim([-20,0])
	plt.xticks(np.arange(-20,0.3,2.5),['-20','','-15','','-10','','-5','','0'])
	plt.margins(0.05,0.05)
	plt.title(sitename, fontsize=12)
	plt.xlabel('Time (kybp)',fontsize=12)
	plt.ylabel('RSL (meters)',fontsize=12)
	#plt.savefig(path.join(ddir_of,'RSL_'+sitename+'.eps'),bbox_inches='tight',pad_inches = 0)
	#pdb.set_trace()
	#plt.show()
	plt.close()

	rsl_all = np.concatenate((site_rsl,rsl_lambeck,rsl_lambeck2),axis=0)
	rsl_time_all = np.concatenate((site_age,age_lambeck,age_lambeck2),axis=0)
	rsl_err_all = np.concatenate((np.stack((np.abs(site_rsl-site_rsl_min),np.abs(site_rsl_max-site_rsl)),axis=0).max(axis=0),rslerr_lambeck,rslerr_lambeck2),axis=0)

	#time for rsl data      
	time_R = -1*rsl_time_all

	#time for case60t
	time_A = topo_rsl01[:,0]
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
	err[i_site] = np.mean(((rsl_all - topo_rsl01[idx_R_60t,1])/rsl_err_all)**2)
	Jacobi[i_site] = np.sin((90.-lat)*np.pi/180.)*1.*np.pi/180.*1.*np.pi/180.
	
#pdb.set_trace()
err_mean = np.sum(Jacobi*err)/np.sum(Jacobi)
print err_mean
