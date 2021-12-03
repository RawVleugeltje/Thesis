# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:55:52 2021

@author: ievdv
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import ERA5_landsea_mask
import netCDF4 as nc

#%% Read data
Wind = xr.open_dataset('Wind/wind.nc')

lats = Wind.latitude.values
lons = Wind.longitude.values

u100 = np.mean(Wind.u100,axis=0)
u10 = np.mean(Wind.u10,axis=0)
v100 = np.mean(Wind.v100,axis=0)
v10 = np.mean(Wind.v10,axis=0)

#%% Wind land sea mask

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')
u100_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, u100)
u10_new = ERA5_landsea_mask.ERA5grid(mask, lons, u10)[0]
v100_new = ERA5_landsea_mask.ERA5grid(mask, lons, v100)[0]
v10_new = ERA5_landsea_mask.ERA5grid(mask, lons, v10)[0]

#%% Calculation of wind direction
wind100 = np.sqrt(u100_new**2 + v100_new**2)
wind10 = np.sqrt(u10_new**2 + v10_new**2)

phi100 = (180 + 180/np.pi * np.arctan2(v100_new,u100_new))%360
phi10 = (180 + 180/np.pi * np.arctan2(v10_new,u10_new))%360

#%% Wind direction to netcdf file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Wind.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

phi100_nc = ncfile.createVariable('Phi100',np.float64,('lat','lon'))
phi100_nc.units = 'degrees'
phi100_nc.standard_name = 'Wind direction at 100 m'

phi10_nc = ncfile.createVariable('Phi10',np.float64,('lat','lon'))
phi10_nc.units = 'degrees' 
phi10_nc.standard_name = 'Wind direction at 10 m'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons_new
phi100_nc[:,:] = phi100
phi10_nc[:,:] = phi10

ncfile.close(); print('Dataset is closed!')