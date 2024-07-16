# %%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import pandas as pd
from scipy.interpolate import griddata
from sklearn.linear_model import LinearRegression

from pyproj import Geod
import cartopy
import cartopy.crs as ccrs

import geoviews as gv
import geoviews.feature as gf
import holoviews as hv
import panel as pn

gv.extension('matplotlib') # Uses MatplotLib backend
gv.output(size=250) # Sets default size for geoview figs
gv.opts.defaults(gv.opts.Layout(sublabel_format="")) # gets rid of subplot lettering for matplotlib combined plots

%matplotlib inline


# %% [markdown]
# ## Initial Data Wrangling 🐴 ##

# %%
# ICE-6G_D Model goes back 122 Ka
ICE6G_D_122Ka = xr.open_dataset('~/VS Code/Zhong Geophysics/ICE_MODELS/ICE6G_122Ka/IceT.I6F_C.131QB_VM5a_1deg.nc')
# ICE-6G_C Model (Updated) Goes back 26 Ka
ICE6G_C_26Ka = xr.open_dataset('~/VS Code/Zhong Geophysics/ICE_MODELS/ICE6G_C_Combined/COMBINED_TIME_ICE-6G_C.nc')

# %%
ICE6G_D_122Ka.load()
print()

# %%
ICE6G_C_26Ka.load()
print()

# %% [markdown]
# ### Constants

# %%
# Approximate m SL/10^6 km^3 ice (meters of sea level rise per 10^6 cubic km of ice)
# From https://linkinghub.elsevier.com/retrieve/pii/S0277379118304074 (sourced from multiple estimates)
mSLRise = [2.485, 2.580, 2.519, 2.478, 2.577, 2.488, 2.466] # m/10^6 km^3

# Volumetric Mean Radius of Earth
R_E = 6371 #km

# Geoid Geometry, using WGS84 Coordinate System (CHECK??)
geod = Geod(ellps="WGS84")


# %% [markdown]
# ### Methods

# %%
 def calc_IceVolume(ice_DS, time_bp):
    """ Estimate Total Volume of ice at a particular time before present.

    Args:
        year (float): thousands of years before present (e.g. input of 26.0 -> 26,000 years before present)
        ice_DS (Dataset): Dataset containing time evolution of ice on Earth
    """
    ice_atyear = ice_DS.sel(Time=time_bp, method="nearest")
    
    total_Volume = 0
    
    ## Need to Find Area of Each Lat/Lon Grid Space
    for lat in ice_atyear["Lat"].values:
        for lon in ice_atyear["Lon"].values:
            # Adjust for lat and lon centering at each 1x1 degree square (e.g. -45.5, 65.5 -> square
            # defined by Lon from -45 to -46 degrees and Lat from 65 to 66 degrees, convert to radians
            lat1 = (lat - 0.5)*np.pi/180 
            lat2 = (lat + 0.5)*np.pi/180
            lon1 = (lon - 0.5)*np.pi/180
            lon2 = (lon + 0.5)*np.pi/180
            
            # Calculate Area of Each square using Archimedes Principle: https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html
            area = (R_E**2)*abs(np.sin(lat1)-np.sin(lat2))*abs(lon1-lon2)
            
            # Find Volume using thickness of ice at this position (stgit in meters) and area
            volume = abs(area)*(ice_atyear.stgit.sel(Lon=lon, Lat=lat).values)/1000 # km^3
            total_Volume += volume
    
    return((total_Volume))

present_day_vol = calc_IceVolume(ICE6G_D_122Ka, 0.0)

# %%
''' CHECK LGM Timing ~26.0 ka
for i in range(20,30):
    print(f"Time: {i} ka -> {round(calc_IceVolume(ICE6G_D_122Ka, i)-present_day_vol,-6)} km^3")

## Expected ~52*10^6 km^3 (Lambeck 2014, https://doi.org/10.1073/pnas.1411762111)
'''

# %% [markdown]
# ## Comparision: Original + Modifications = Proposed

# %% [markdown]
# ### Plotting Comparison Methods

# %%
def plot_ice_dist(model, time):
    ice_test = model.stgit # thickness variable

    fig = plt.subplots(figsize=(9,6))

    ax = plt.axes(projection = ccrs.Mercator())
    ax.coastlines()
    ax.gridlines()
    
    
    plot_elem = ice_test.sel(Time = time, method = "nearest").plot(x = 'Lon', y = 'Lat', transform = ccrs.PlateCarree(), cmap="RdBu_r", cbar_kwargs={'shrink':0.4}, ax = ax)
    

# %%
def plot_ice_comparison(model_orig, model_new, time, new_label='Proposed'):
    """Plot the three datasets side by side allowing us to see the original distribution, changes we made, and the final map

    Args:
        model_orig (DataArray): Original Distribution of Ice
        model_new (DataArray): New modified distribution of ice
        time (float): number of thousands of years ago that we would like to access (e.g. Input:26.0 -> 26,000 years ago)
        proposed_label (string): identifier for proposed model
    """
    model_difference = (model_new - model_orig)*20
    
    model_triple = xr.concat([model_orig, model_difference, model_new], pd.Index(['Original', 'Difference (scaled by 10)', new_label],name='Version'))

    
    plot_elem = model_triple.stgit.sel(Time = time, method = "nearest").isel(Version = slice(0,3)).plot(x = 'Lon', y = 'Lat',col = 'Version',transform = ccrs.PlateCarree(), cmap="RdBu_r", cbar_kwargs={'shrink':0.4}, subplot_kws={'projection':ccrs.Mercator()}, aspect = 1, size = 7)
    
    for ax in plot_elem.axs.flat:
        ax.coastlines()
        ax.gridlines()
    
    
    return model_triple
        
        
    

# %%
def plot_ice_comparison_alt(model_orig, model_new, time, new_label='Proposed'):
    """Plot the three datasets side by side allowing us to see the original distribution, changes we made, and the final map
       Changes from original plot_ice_comparison() -> alt version allows for custom scaling and coloring of each graph and other attributes

    Args:
        model_orig (DataArray): Original Distribution of Ice
        model_new (DataArray): New modified distribution of ice
        time (float): number of thousands of years ago that we would like to access (e.g. Input:26.0 -> 26,000 years ago)
        proposed_label (string): identifier for proposed model
    """
    ## Change in ice thickness between models
    model_diff = model_new - model_orig
    
    ## Convert datasets to geoview datasets
    gv_orig = gv.Dataset(model_orig.copy().sel(Time = 26.0), ['Lon', 'Lat'], 'stgit', crs=ccrs.PlateCarree())
    gv_new = gv.Dataset(model_new.copy().sel(Time = 26.0), ['Lon', 'Lat'], 'stgit', crs=ccrs.PlateCarree())
    gv_diff = gv.Dataset(model_diff.copy().sel(Time = 26.0), ['Lon', 'Lat'], 'stgit', crs=ccrs.PlateCarree())
    
    ## Generate geoview images for each dataset
    image_orig = gv_orig.to(gv.Image)
    image_new = gv_new.to(gv.Image)
    image_diff = gv_diff.to(gv.Image)
    
    ## Calibrate options for each image and plot
    plot_orig = image_orig.opts(projection=ccrs.Mercator(),cmap='RdBu_r', colorbar=True, title = 'Original', xticks=5, yticks = 7)*gf.coastline
    plot_new = image_new.opts(projection=ccrs.Mercator(),cmap='RdBu_r', colorbar=True, title = new_label, xticks=5, yticks = 7)*gf.coastline
    plot_diff = image_diff.opts(projection=ccrs.Mercator(),cmap='Viridis', colorbar=True, clim = (model_diff.stgit.min(), model_diff.stgit.max()), title = 'Difference', xticks=5, yticks = 7)*gf.coastline ## Note: The options for this one suck and are very long, mostly so that the color bar has custom scaling, see clim attr
    
    
    return plot_orig+plot_diff+plot_new
        
        
    

# %% [markdown]
# ## Different Proposed Models

# %% [markdown]
# #### Original Model

# %%
time = 26.0 # Ka

print("Original Ice Volume: ", round(calc_IceVolume(ICE6G_D_122Ka,time) - present_day_vol,-6))

# %% [markdown]
# #### Simple Linear Model

# %%
# Modification Constants
m_const = 1.01
b_const = 0

Linear_Test = ICE6G_D_122Ka.copy()
Linear_Test['stgit'] = Linear_Test.stgit*m_const + b_const*(Linear_Test.stgit > 200)

print("Simple Linear Ice Volume: ", round(calc_IceVolume(Linear_Test, time) - present_day_vol,-6))

# %% [markdown]
# #### Latitude-Dependent Linear Model

# %%
# Modifying Constants
m_lat_north = 1.04
b_lat_north = 0
m_lat_south = 0.98
b_lat_south = 0

Linear_Lat_Test = ICE6G_D_122Ka.copy()

# Retrieve Ice thickness separately
northern_stgit = Linear_Lat_Test['stgit'].sel(Lat=slice(0, 90)) # Including equator in northern and not southern, shouldn't matter at all though
southern_stgit = Linear_Lat_Test['stgit'].sel(Lat=slice(-90,0))

# Modify with set constants
northern_stgit = northern_stgit*m_lat_north + b_lat_north
southern_stgit = southern_stgit*m_lat_south + b_lat_south

Linear_Lat_Test = xr.merge([northern_stgit, southern_stgit])

print("Latitude-Dependent Linear Ice Volume: ", round(calc_IceVolume(Linear_Lat_Test, time) - present_day_vol,-6))

# %% [markdown]
# #### Quadratic Model

# %%
# Modifying Constants
second_deg = 0.000043
first_deg = -0.1
zero_deg = 0

Quadratic_Test = ICE6G_D_122Ka.copy()

# Retrieve Ice thickness
Quad_stgit = Quadratic_Test['stgit']

# Modify with set constants
Quad_stgit = Quad_stgit + (Quad_stgit**2)*second_deg + (Quad_stgit)*first_deg + zero_deg

Quadratic_Test['stgit'] = Quad_stgit

print("Latitude-Dependent Quadratic Ice Volume: ", round(calc_IceVolume(Quadratic_Test, time) - present_day_vol,-6))

# %% [markdown]
# #### Latitude-Dependent Quadratic Model

# %%
# Modifying Constants
second_lat_north = 0.000043
first_lat_north = -0.06
zero_lat_north = 0
second_lat_south = -0.00002
first_lat_south = 0
zero_lat_south = 0

Quadratic_Lat_Test = ICE6G_D_122Ka.copy()

# Retrieve Ice thickness separately
northern_stgit = Quadratic_Lat_Test['stgit'].sel(Lat=slice(0, 90)) # Including equator in northern and not southern, shouldn't matter at all though
southern_stgit = Quadratic_Lat_Test['stgit'].sel(Lat=slice(-90,0))

# Modify with set constants
northern_stgit = northern_stgit + (northern_stgit**2)*second_lat_north + northern_stgit*first_lat_north + zero_lat_north
southern_stgit = southern_stgit + (southern_stgit**2)*second_lat_south + southern_stgit*first_lat_south + zero_lat_south

Quadratic_Lat_Test = xr.merge([northern_stgit, southern_stgit])

print("Latitude-Dependent Quadratic Ice Volume: ", round(calc_IceVolume(Quadratic_Lat_Test, time) - present_day_vol,-6))

# %% [markdown]
# #### Future Methods ####
# - Only NA Ice Sheets
# - Latitude-Dependent Gradient
# - Focus even more around hudson bay
# - Superimpose ANU Model (https://www.pnas.org/doi/abs/10.1073/pnas.1411762111)
# - Look at replacing NA ice sheets with Gowan 2021 (https://www.nature.com/articles/s41467-021-21469-w)
# 

# %% [markdown]
# ## Visualizing Different Proposed Methods

# %% [markdown]
# #### Simple Linear Model

# %%
plot_ice_comparison_alt(ICE6G_D_122Ka, Linear_Test, 26.0, 'Simple Linear Model')

# %% [markdown]
# #### Latitude-Dependent Linear Model

# %%
plot_ice_comparison_alt(ICE6G_D_122Ka, Linear_Lat_Test, 26.0, 'Latitude-Dependent Linear Model')

# %% [markdown]
# #### Quadratic Model

# %%
plot_ice_comparison_alt(ICE6G_D_122Ka, Quadratic_Test, 26.0, 'Quadratic Model')

# %% [markdown]
# #### Latitude-Dependent Quadratic Model

# %%
plot_ice_comparison_alt(ICE6G_D_122Ka, Quadratic_Lat_Test, 26.0, 'Latitude-Dependent Quadratic Model')

# %%



