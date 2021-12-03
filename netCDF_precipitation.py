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

#%% Read data
PrecTemp = xr.open_dataset('Prec_Temp/prectemp.nc')

lats = PrecTemp.latitude.values
lons = PrecTemp.longitude.values

a = np.sum(PrecTemp.tp[0:12,:,:],axis=0)*1000

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')
a_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, a)

plt.title('Methode 1')
plt.imshow(a_new)
plt.colorbar()
plt.show()

fig, ax = plt.subplots(1,2,figsize=(15,5))
ax[0].set_title('Methode 1 volledige colorbar')
a2plot1 = ax[0].imshow(a_new)
plt.colorbar(a2plot1, ax=ax[0], orientation='horizontal')
ax[1].set_title('Methode 1 gelimiteerde colorbar')
a2plot2 = ax[1].imshow(a_new,vmax=100)
plt.colorbar(a2plot2, ax=ax[1], orientation='horizontal')
plt.show()

days = [31,28,31,30,31,30,31,31,30,31,30,31]
a2 = PrecTemp.tp[0:12,:,:]
for i in range(12):
    a2[i] = PrecTemp.tp[i,:,:]*days[i]*1000

a2 = np.sum(a2[0:12,:,:],axis=0)
a2_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, a2)

fig, ax = plt.subplots(1,2,figsize=(15,5))
ax[0].set_title('Methode 2 volledige colorbar')
a2plot1 = ax[0].imshow(a2_new)
plt.colorbar(a2plot1, ax=ax[0], orientation='horizontal')
ax[1].set_title('Methode 2 gelimiteerde colorbar')
a2plot2 = ax[1].imshow(a2_new,vmax=3000)
plt.colorbar(a2plot2, ax=ax[1], orientation='horizontal')
plt.show()
#%% Precipitation land sea mask

# #%% Precipitation netCDF file
# try: ncfile.close()
# except: pass
# ncfile = nc.Dataset('Precipitation.nc',mode='w',format='NETCDF4_CLASSIC') 
# print(ncfile)

# lat_dim = ncfile.createDimension('lat', len(lats))
# lon_dim = ncfile.createDimension('lon', len(lons))
# month_dim = ncfile.createDimension('month', 12)

# lat = ncfile.createVariable('lat', np.float32, ('lat',))
# lat.units = 'degrees'
# lat.long_name = 'latitude'
# lon = ncfile.createVariable('lon', np.float32, ('lon',))
# lon.units = 'degrees'
# lon.long_name = 'longitude'
# month = ncfile.createVariable('month', np.float32, ('month'))
# month.units = '-'
# month.long_name = 'month'

# temp_mean = ncfile.createVariable('T_mean',np.float64,('lat','lon'))
# temp_mean.units = 'K'
# temp_mean.standard_name = 'Mean temperature'

# temp_monthly = ncfile.createVariable('T_monthly',np.float64,('month','lat','lon'))
# temp_monthly.units = 'K'
# temp_monthly.standard_name = 'Monthly mean temperature'

# nlats = len(lat_dim); nlons = len(lon_dim)
# lat[:] = lats
# lon[:] = lons
# month[:] = [1,2,3,4,5,6,7,8,9,10,11,12]
# temp_mean[:,:] = t2m_new
# temp_monthly[:,:,:] = t2m_monthly_new

# ncfile.close(); print('Dataset is closed!')