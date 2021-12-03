# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:23:15 2021

@author: ievdv
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import ERA5_landsea_mask
import netCDF4 as nc

#%% Read data
PrecTemp = xr.open_dataset('Prec_Temp/prectemp.nc')

lats = PrecTemp.latitude.values
lons = PrecTemp.longitude.values

#%% Temperature land sea mask

t2m = np.mean(PrecTemp['t2m'],axis=0)
t2m_monthly = PrecTemp['t2m'].groupby('time.month').mean()

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')
t2m_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, t2m)

t2m_monthly_new = t2m_monthly
for i in range(12):
    t2m_monthly_new[i,:,:] = ERA5_landsea_mask.ERA5grid(mask, lons, t2m_monthly[i,:,:])[0]
    
#%% Temperature netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Temperature.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))
month_dim = ncfile.createDimension('month', 12)

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'
month = ncfile.createVariable('month', np.float32, ('month'))
month.units = '-'
month.long_name = 'month'

temp_mean = ncfile.createVariable('T_mean',np.float64,('lat','lon'))
temp_mean.units = 'K'
temp_mean.standard_name = 'Mean temperature'

temp_monthly = ncfile.createVariable('T_monthly',np.float64,('month','lat','lon'))
temp_monthly.units = 'K'
temp_monthly.standard_name = 'Monthly mean temperature'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons_new
month[:] = [1,2,3,4,5,6,7,8,9,10,11,12]
temp_mean[:,:] = t2m_new
temp_monthly[:,:,:] = t2m_monthly_new

ncfile.close(); print('Dataset is closed!')