#!/usr/bin/env python 
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas
from scipy.interpolate import interp2d

#pdb.set_trace()
rmsData = pandas.read_excel('/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/gridresearch_results/RMS_viscosity_grids.xlsx')
viscUpperMantle = rmsData['upper mantle visc'].values
viscLowerMantle = rmsData['lower mantle visc'].values
misfit = rmsData['misfit_ANU_Lambeck_na'].values
#misfit = rmsData['misfit_N_F_BAR'].values
#misfit = rmsData['misfit_lambeck_na'].values
#misfit = rmsData['misfit_lambeck_f'].values
#misfit = rmsData['misfit_N_F'].values
#misfit = rmsData['misfit_peltier'].values

plt.figure()
#pdb.set_trace()
#plt.pcolor(viscLowerMantle.reshape([6, 7]), viscUpperMantle.reshape([6, 7]), misfit.reshape([6, 7]), cmap='jet')
x0_grid = viscLowerMantle.reshape([6, 7])
y0_grid = viscUpperMantle.reshape([6, 7])
z0_grid = misfit.reshape([6, 7]) 
f = interp2d(x0_grid, y0_grid, z0_grid, kind='linear')

x = np.arange(20.5, 23.6, 0.1)
y = np.arange(19, 21.6, 0.1)
x_grid = np.array([x] * len(y))
y_grid = np.array([y] * len(x))
y_grid = y_grid.T
z_grid = f(x, y)


#plt.contourf(x_grid, y_grid, z_grid, cmap='jet')
#plt.pcolor(x_grid, y_grid, z_grid, cmap='jet')

x00 = np.arange(20.5-0.25, 23.6+0.25, 0.5)
y00 = np.arange(19-0.25, 21.6+0.25, 0.5)
x00_grid = np.array([x00] * len(y00))
y00_grid = np.array([y00] * len(x00))
y00_grid = y00_grid.T
plt.pcolormesh(x00_grid, y00_grid, z0_grid, shading='flat', cmap='jet')
plt.xlabel('Lower mantle viscosity (log)')
plt.ylabel('Upper mantle viscosity (log)')
#plt.yticks(np.arange(19.0, 22, 0.5))
plt.colorbar()
#plt.title('Peltier global dataset(954 sites)')
#plt.title('Peltier 2015 paper dataset(30 sites)')
#plt.title('Lambeck far field dataset(36 sites)')
plt.title('ANU VS Lambeck North America dataset(64 sites)')
plt.savefig('rms_viscosity_ANU_Lambeck_NA.png')
plt.show()

