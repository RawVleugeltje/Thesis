# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:12:55 2021

@author: ievdv
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import ERA5_landsea_mask
import netCDF4 as nc
import glob

#%% Read data
files_prec = glob.glob('Precipitation/era5_total_precipitation_*.nc')

prec = []
for i in range(len(files_prec)):
    prec.append(xr.open_dataset(files_prec[i]).tp[0,:,:])

prec_mean = np.mean(prec,axis=0)
lats = xr.open_dataset(files_prec[0]).latitude
lons = xr.open_dataset(files_prec[0]).longitude

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')
prec_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, prec_mean)

#%% Precipitation netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Precipitation.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

prec_mean = ncfile.createVariable('Prec_mean',np.float64,('lat','lon'))
prec_mean.units = 'm'
prec_mean.standard_name = 'Mean precipitation'

lat[:] = lats
lon[:] = lons_new
prec_mean[:,:] = prec_new

ncfile.close(); print('Dataset is closed!')