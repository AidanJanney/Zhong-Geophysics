#!/usr/bin/env python
u"""
ncdf_write.py
Written by Tyler Sutterley (UPDATED: 07/2016)

This program writes COARDS-compliant NetCDF files for GMT and other programs

CALLING SEQUENCE:
	ncdf_write(data, lon, lat, time, FILENAME='sigma.nc', TITLE = 'Spatial Data',
		LONGNAME='Equivalent Water Thickness', UNITS='cmH2O/yr')
INPUTS:
	data: z data
	lon: longitude array
	lat: latitude array
	time: date of z data
OPTIONS:
	FILENAME: output filename netcdf
	UNITS: z variable units (cmH2O, etc.)
	LONGNAME: z variable description (Equivalent Water Thickness, etc.)
	VARNAME: z variable name in netcdf file
	LONNAME: longitude variable name in netcdf file
	LATNAME: latitude variable name in netcdf file
	TIME_UNITS: time variable units (default: years)
	TIME_LONGNAME: time variable description (default: time)
	FILL_VALUE: missing value for z variable
	TITLE: title attribute of dataset
	CLOBBER: will overwrite an existing netcdf file
	VERBOSE: will print to screen the netcdf structure parameters

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)
	netCDF4: Python interface to the netCDF C library
	 	(http://unidata.github.io/netcdf4-python/)

UPDATE HISTORY:
	Updated 10/2016: formating input to be data[time,lat,lon] for Meng Zhao's tradition
	Updated 07/2016: using netCDF4-python with zlib compression
	Updated 06/2016: using __future__ print function
	Updated 05/2016: output data types same as input data types
	Updated 11/2014: new parameters for variable names and attributes
	Updated 05/2014: new parameters for time attributes, and missing values
	Updated 02/2014: minor update to if statements
	Updated 07/2013: switched from Scientific Python to Scipy
	Updated 04/2013: converted to python
	Updated 03/2013: converted to Octave
	Updated 01/2013: adding time as a variable
	Updated 10/2012: changed from variable names x and y to lon and lat.
		As xyz2grd used x and y, I thought GMT needed x and y.
		Other programs also use lon-lat.  This makes more sense anyways
	Written 07/2012 for GMT and for archiving datasets
		Motivation for archival: NetCDF files are much smaller than ascii
		files and more portable/transferable than IDL .sav files
		(possible to connect with geostatistics packages in R?)
"""
from __future__ import print_function

import netCDF4
import numpy as np

def ncdf_write(data, lon, lat, time, FILENAME='sigma.nc', UNITS='cmH2O', \
	LONGNAME='Equivalent Water Thickness', LONNAME='lon', LATNAME='lat', \
	VARNAME='z', TIME_UNITS='years', TIME_LONGNAME='time', FILL_VALUE=None, \
	TITLE = 'Spatial Data', CLOBBER='Y', VERBOSE='N'):

	#-- setting NetCDF clobber attribute
	if CLOBBER in ('Y','y'):
		clobber = 'w'
	else:
		clobber = 'a'

	#-- opening NetCDF file for writing
	#-- Create the NetCDF file
	fileID = netCDF4.Dataset(FILENAME, clobber, format="NETCDF4")

	#-- Dimensions of parameters
	n_lon = len(lon)
	n_lat = len(lat)
	if (np.ndim(time) == 0):
		n_time = 1
	else:
		n_time = len(time)

	#-- Defining the NetCDF dimensions
	fileID.createDimension(LONNAME, n_lon)
	fileID.createDimension(LATNAME, n_lat)
	fileID.createDimension('time', n_time)

	#-- defining the NetCDF variables
	#-- lat and lon
	nc_lon = fileID.createVariable(LONNAME, lon.dtype, (LONNAME,))
	nc_lat = fileID.createVariable(LATNAME, lat.dtype, (LATNAME,))
	#-- spatial data
	if (n_time > 1):
		nc_z = fileID.createVariable(VARNAME, data.dtype,
			('time',LATNAME,LONNAME,), fill_value=FILL_VALUE, zlib=True)
	else:
		nc_z = fileID.createVariable(VARNAME, data.dtype,
			(LATNAME,LONNAME,), fill_value=FILL_VALUE, zlib=True)
	#-- time (in decimal form)
	nc_time = fileID.createVariable('time', 'f8', ('time',))

	#-- filling NetCDF variables
	nc_lon[:] = lon
	nc_lat[:] = lat
	if (n_time > 1):
		nc_z[:,:,:] = data
	else:
		nc_z[:,:] = data
	nc_time[:] = time

	#-- Defining attributes for longitude and latitude
	nc_lon.long_name = 'longitude'
	nc_lon.units = 'degrees_east'
	nc_lat.long_name = 'latitude'
	nc_lat.units = 'degrees_north'
	#-- Defining attributes for dataset
	nc_z.long_name = LONGNAME
	nc_z.units = UNITS
	#-- Defining attributes for date
	nc_time.long_name = TIME_UNITS
	nc_time.units = TIME_LONGNAME
	#-- global variable of NetCDF file
	fileID.TITLE = TITLE

	#-- Output NetCDF structure information
	if VERBOSE in ('Y','y'):
		print(FILENAME)
		print(fileID.variables.keys())

	#-- Closing the NetCDF file
	fileID.close()
