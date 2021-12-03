# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 09:00:23 2021

@author: ievdv
"""
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import LandSea_mask
import numpy as np
import glob
import netCDF4 as nc

#%%
era5 = xr.open_dataset('ERA5/ERA5-HRES_MEV_results.nc')
gpex = xr.open_dataset('GPEX/GPEX.nc')
era5_T1000 = np.loadtxt('ERA5/T1000.csv')
era5_T100 = np.loadtxt('ERA5/T100.csv')
era5_T10 = np.loadtxt('ERA5/T10.csv')
era5_T1 = np.loadtxt('ERA5/T1.csv')
era5_T5 = np.loadtxt('ERA5/T5.csv')

#%%
lats = era5.lat.values
lons = era5.lon.values

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')

var10 = era5['quants'][:,:,0]
var100 = era5['quants'][:,:,3]
var50 = era5['quants'][:,:,2]
var500 = era5['quants'][:,:,4]

rlon = len(lons)
rlonhalf = int(rlon/2)

# east west verkeerd
maskeast = mask[:,rlonhalf:]
maskwest = mask[:,:rlonhalf]
maskr    = np.hstack((maskeast, maskwest))

varmask10 = var10 * maskr
varmask100 = var100 * maskr
varmask50 = var50 * maskr
varmask500 = var500 * maskr

lons_new = lons[lons>=180]-360
lons_old = lons[lons<180]
lons = np.concatenate((lons_old,lons_new))

var_new10 = np.zeros_like(varmask10)
var_new100 = np.zeros_like(varmask100)
var_new50 = np.zeros_like(varmask50)
var_new500 = np.zeros_like(varmask500)
lons_new = np.zeros_like(lons)

for i in range(len(lons)):
    if lons[i] < 0:
        var_new10[:,i-720] = varmask10[:,i]
        var_new100[:,i-720] = varmask100[:,i]
        var_new50[:,i-720] = varmask50[:,i]
        var_new500[:,i-720] = varmask500[:,i]
        lons_new[i-720] = lons[i]
    else:
        var_new10[:,i+720] = varmask10[:,i]
        var_new100[:,i+720] = varmask100[:,i]
        var_new50[:,i+720] = varmask50[:,i]
        var_new500[:,i+720] = varmask500[:,i]
        lons_new[i+720] = lons[i]
        
# # plot T5
# fig, ax = plt.subplots(1,2,figsize=(12,10))
# var1 = ax[0].imshow(gpex['mev_estimate'][:,:,3,1],vmax=100)
# ax[0].set_title('T5 GPEX')
# plt.colorbar(var1,ax=ax[0],orientation='horizontal',pad=0.05)

# var2 = ax[1].imshow(era5_T5[:610,:],vmax=100)
# ax[1].set_title('T5 calculated from era5')
# plt.colorbar(var1,ax=ax[1],orientation='horizontal',pad=0.05)


# # plot T10
# fig, ax = plt.subplots(3,1,figsize=(12,10))
# var1 = ax[0].imshow(var_new,extent=[lons.min(),lons.max(),lats.min(),
#                                     lats.max()],vmax=200)
# ax[0].set_title('T10 era5')
# # plt.colorbar(var1,ax=ax[0],orientation='horizontal',pad=0.1,fraction=0.1)

# var2 = ax[1].imshow(era5_T10,vmax=200)
# ax[1].set_title('T10 calculated from era5')
# # plt.colorbar(var1,ax=ax[1],orientation='horizontal',pad=0.1,fraction=0.1)

# var3 = ax[2].imshow(gpex['mev_estimate'][:,:,3,2],vmax=200)
# ax[2].set_title('T10 GPEX')
# plt.colorbar(var1,ax=ax[2],orientation='horizontal',pad=0.1,fraction=0.1)
# plt.tight_layout()
# plt.show()

# # plot T100
# fig, ax = plt.subplots(1,2,figsize=(12,10))
# var1 = ax[0].imshow(var_new100,extent=[lons.min(),lons.max(),lats.min(),
#                                     lats.max()],vmax=500)
# ax[0].set_title('T100 era5')
# plt.colorbar(var1,ax=ax[0],orientation='horizontal',pad=0.05)

# var2 = ax[1].imshow(era5_T100,vmax=500)
# ax[1].set_title('T100 calculated from era5')
# plt.colorbar(var1,ax=ax[1],orientation='horizontal',pad=0.05)
# plt.show()

# # plot T1000
# fig, ax = plt.subplots(1,2,figsize=(12,10))
# var1 = ax[0].imshow(gpex['mev_estimate'][:,:,3,9],vmax=1000)
# ax[0].set_title('T1000 GPEX')
# plt.colorbar(var1,ax=ax[0],orientation='horizontal',pad=0.05)

# var2 = ax[1].imshow(era5_T1000[:610,:],vmax=1000)
# ax[1].set_title('T1000 calculated from era5')
# plt.colorbar(var1,ax=ax[1],orientation='horizontal',pad=0.05)

#%%
h5 = (var_new500 - 2*var_new50 + era5_T5) / (var_new50 - era5_T5)
h1 = (var_new100 - 2*var_new10 + era5_T1) / (var_new10 - era5_T1)
h10 = (era5_T1000 - 2*var_new100 + var_new10) / (var_new100 - var_new10)

h5 = h5
h1 = h1
h10 = h10

lats = era5.lat.values
lons = lons_new

fig, ax = plt.subplots(2,2,figsize=(16,8))
plot1 = ax[0,0].imshow(h1,vmin=-1.0,vmax=1.0)
plt.colorbar(plot1,ax=ax[0,0],extend='both',orientation='horizontal')
ax[0,0].set_title('Heaviness calculated with T100,T10,T1 from era5')

plot2 = ax[0,1].imshow(h5,vmin=-1.0,vmax=1.0)
plt.colorbar(plot2,ax=ax[0,1],extend='both',orientation='horizontal')
ax[0,1].set_title('Heaviness calculated with T500,T50,T5 from era5')

plot4 = ax[1,0].imshow(h10,vmin=-1.0,vmax=1.0)
plt.colorbar(plot2,ax=ax[1,0],extend='both',orientation='horizontal')
ax[1,0].set_title('Heaviness calculated with T1000,T100,T10 from era5')


plot3 = ax[1,1].imshow(gpex['mev_heaviness'][:,:,3],vmin=-1.0,vmax=1.0)
plt.colorbar(plot3,ax=ax[1,1],extend='both',orientation='horizontal')
ax[1,1].set_title('Heaviness calculated with T1000,T100,T10 from GPEX')
plt.tight_layout()
plt.show()

#%% make netcdf4
try: ncfile.close()
except: pass
ncfile = nc.Dataset('Heaviness_H1_H5_H10.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(lats))
lon_dim = ncfile.createDimension('lon', len(lons))

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'
heav1 = ncfile.createVariable('h1',np.float64,('lat','lon'))
heav1.units = '-'
heav1.standard_name = 'heaviness T1, T10, T100'

heav5 = ncfile.createVariable('h5',np.float64,('lat','lon'))
heav5.units = '-'
heav5.standard_name = 'heaviness T5, T50, T500'

heav10 = ncfile.createVariable('h10',np.float64,('lat','lon'))
heav10.units = '-'
heav10.standard_name = 'heaviness T10, T100, T1000'

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = lats
lon[:] = lons
heav1[:,:] = h1
heav5[:,:] = h5
heav10[:,:] = h10


ncfile.close(); print('Dataset is closed!')