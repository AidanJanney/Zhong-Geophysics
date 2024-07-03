#/usr/bin/env python
#this code is used to generate the (lon lat) file for the ice covered region at LGM

from os import path
import numpy as np
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pdb

ddir = '/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/ice6g_data'
fn = 'I6_C.VM5a_1deg.26.nc'
data = netCDF4.Dataset(path.join(ddir,fn),'r')
#lon = np.array(data['lon'])
#lat = np.array(data['lat'])
#ih = np.array(data['stgit'])

#pdb.set_trace()
#lon = lon-180.0
#ih1 = ih

#plot subregion in North American
data1=np.loadtxt('Ameri_subset.txt')
lon = data1[:,0]
lat = data1[:,1]
ih = data1[:,2]

plt.figure()
plt.plot(lon,lat,'.')
plt.show()

#pdb.set_trace()
lon = lon -180.0
shp = np.shape(ih)
ih1 = np.zeros((shp[0],shp[1]))
for i_lon in range(len(lon)):
    if i_lon <= 179:
        ih1[:,i_lon] = ih[:,i_lon+179]
for i_lon in range(len(lon)):
    if i_lon > 179:
        ih1[:,i_lon] = ih[:,i_lon-180]

#plote map of uplift global distribution
lon_0 = lon.mean()
lat_0 = lat.mean()
#m = Basemap(width=5000000,height=3500000,resolution='l',projection='stere',\
#            lat_0=lat_0,lon_0=lon_0)
m = Basemap()
#lons,lats = np.meshgrid(lon,lat)
#xi,yi = m(lons,lats)

#cs = m.pcolor(xi,yi,ih1)
cs = m.pcolor(lon,lat,ih)
#m.drawparallels(np.arange(-90.,90.,20.))
#m.drawmeridians(np.arange(-180,180,30))
m.drawcoastlines()
#plt.colorbar()
plt.style.use('classic')
plt.colorbar(cs,orientation='horizontal',label='meter')
#plt.title('surface geoid rate at present day')
#plt.savefig('map_geoid_rate_presentday')
plt.title('ice distribution at LGM')
#plt.savefig('ice_map_LGM')
plt.show()

