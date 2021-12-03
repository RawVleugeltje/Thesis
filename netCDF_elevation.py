# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 14:21:08 2021

@author: ievdv
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import ERA5_landsea_mask
import netCDF4 as nc
import pandas as pd
from scipy.interpolate import griddata

#%% Read data
elevation = xr.open_dataset('Elevation/GMTED2010_15n060_0250deg.nc')

lats = elevation.latitude.values
lons = elevation.longitude.values

#%% Elevation land sea mask

elev = elevation.elevation

mask = LandSea_mask.create_mask_array(lats,lons,shapefile='continent',ysign=1)

elev_new = elev*np.flipud(mask)

#%% Interpolation to ERA5 grid
era5 = xr.open_dataset('ERA5/ERA5-HRES_MEV_results.nc')

lat_era5 = era5.lat.values
lon_era5 = era5.lon.values - 180

grid_x = np.full((len(lon_era5),len(lat_era5)),lat_era5)
grid_y = np.full((len(lat_era5),len(lon_era5)),lon_era5)

grid_x = np.transpose(grid_x)

values = elev_new.values.flatten()
points = np.zeros((len(values),2))

a = []
for i in range(len(lats)):
    for j in range(len(lons)):
        a.append(i*1440+j)
        points[i*1440+j] = lats[i], lons[j]

elev_era5 = griddata(points,values,(grid_x,grid_y),method='linear')

#%% Elevation netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Elevation.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lat_era5))
lon_dim = ncfile.createDimension('lon', len(lon_era5))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

elevation = ncfile.createVariable('Elevation',np.float64,('lat','lon'))
elevation.units = 'm'
elevation.standard_name = 'Mean elevation in m'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lat_era5
lon[:] = lon_era5
elevation[:,:] = elev_era5

ncfile.close(); print('Dataset is closed!')