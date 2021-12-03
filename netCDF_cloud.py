# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:29:41 2021

@author: ievdv
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import ERA5_landsea_mask
import netCDF4 as nc

#%% Read data
cloud = xr.open_dataset('CloudFraction/cloudfrac19792019.nc')

hcc = np.mean(cloud.hcc[:,0,:,:],axis=0)
mcc = np.mean(cloud.mcc[:,0,:,:],axis=0)
lcc = np.mean(cloud.lcc[:,0,:,:],axis=0)
tcc = np.mean(cloud.tcc[:,0,:,:],axis=0)

lats = cloud.latitude.values
lons = cloud.longitude.values

#%% Cloud land sea mask

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')

hcc_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, hcc)
mcc_new = ERA5_landsea_mask.ERA5grid(mask, lons, mcc)[0]
lcc_new = ERA5_landsea_mask.ERA5grid(mask, lons, lcc)[0]
tcc_new = ERA5_landsea_mask.ERA5grid(mask, lons, tcc)[0]

#%% Temperature netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Clouds.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

hcc = ncfile.createVariable('hcc',np.float64,('lat','lon'))
hcc.units = '-'
hcc.standard_name = 'High cloud cover'

mcc = ncfile.createVariable('mcc',np.float64,('lat','lon'))
mcc.units = '-'
mcc.standard_name = 'Middle cloud cover'

lcc = ncfile.createVariable('lcc',np.float64,('lat','lon'))
lcc.units = '-'
lcc.standard_name = 'Low cloud cover'

tcc = ncfile.createVariable('tcc',np.float64,('lat','lon'))
tcc.units = '-'
tcc.standard_name = 'Total cloud cover'


nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons_new
hcc[:,:] = hcc_new
mcc[:,:] = mcc_new
lcc[:,:] = lcc_new
tcc[:,:] = tcc_new

ncfile.close(); print('Dataset is closed!')