# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:44:50 2021

@author: ievdv
"""

import xarray as xr
import rioxarray as rio
import richdem as rd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

#%% make tif file
elev = xr.open_dataset('netCDF/Elevation.nc')

lats = elev.lat
lons = elev.lon

elev = elev.rename({'lat':'y'})
elev = elev.rename({'lon':'x'})
elev.rio.to_raster('elevation.tif')

#%% open tiffile and calculate aspect
dem = rd.LoadGDAL('elevation.tif', no_data=-9999)

aspect = rd.TerrainAttribute(dem, attrib='aspect')
rd.rdShow(aspect, cmap='viridis', figsize=(8, 5.5))
plt.show()

#%% Aspect netCDF file
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Slope.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'

aspect_nc = ncfile.createVariable('Aspect',np.float64,('lat','lon'))
aspect_nc.units = 'degrees'
aspect_nc.standard_name = 'Aspect of the slope'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons
aspect_nc[:,:] = aspect

ncfile.close(); print('Dataset is closed!')