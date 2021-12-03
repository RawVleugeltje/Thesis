# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:34:22 2021

@author: ievdv
"""
import xarray as xr
import numpy as np
import netCDF4 as nc

#%% Import data
slope = xr.open_dataset('netCDF/Slope.nc')
wind = xr.open_dataset('netCDF/Wind.nc')

lats = slope.lat
lons = slope.lon

#%% Calculate slope/wind proxy
rho_int = np.abs(wind.Phi100.values - slope.Aspect.values)

rho = rho_int

for i in range(rho.shape[0]):
    print(i/721*100)
    for j in range(rho.shape[1]):
        if rho_int[i,j] <= 180:
            rho[i,j] = 1 - rho_int[i,j] / 180
        else:
            rho[i,j] = 1 - (360 - rho_int[i,j]) / 180
            
#%%
#%% Wind slope to netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('windslope.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

windslope = ncfile.createVariable('WindSlope',np.float64,('lat','lon'))
windslope.units = '-'
windslope.standard_name = 'rho, representing the wind/slope'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons
windslope[:,:] = rho

ncfile.close(); print('Dataset is closed!')
