# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 11:43:34 2021

@author: ggrundemann
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import regionmask
import netCDF4 as nc
import matplotlib

#%% open dataset
h = xr.open_dataset('netCDF/Heaviness_H1_H5_H10.nc')
GPEX = xr.open_dataset('heaviness/GPEX.nc')

hGPEX = GPEX.mev_heaviness[:,:,3]
#%% all IPCC reference regions 
ar6_all = regionmask.defined_regions.ar6.all
ar6_all

mask_3D = ar6_all.mask_3D(h)
mask_3D

mask_3D_GPEX = ar6_all.mask_3D(hGPEX)
mask_3D_GPEX

#%% land mask
land_110 = regionmask.defined_regions.natural_earth.land_110

land_mask = land_110.mask_3D(h)
land_mask_GPEX = land_110.mask_3D(hGPEX)

#%% quick plot
ax = plt.axes(projection=ccrs.Robinson())
land_mask.sel(region=0).plot.imshow(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
plt.show()

ax = plt.axes(projection=ccrs.Robinson())
land_mask_GPEX.sel(region=0).plot.imshow(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
plt.show()

#%% plot of all the regions seperatly
plot_new = np.zeros_like(h.h10)

for i in range(len(mask_3D)):
    a = mask_3D.isel(region=i).values
    for q in range(721):
        for p in range(1440):
            if a[q,p] == True:
                plot_new[q,p] = 1*i
                
plot_lat = land_mask.lat.values
plot_lon = land_mask.lon.values

#%%
from collections import defaultdict

a_color = defaultdict(list)
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
for i in range(58):
     a_color[str(i)].append(np.array(matplotlib.colors.to_rgb(cmaplist[4*i])))
     
data_3d = np.ndarray(shape=(plot_new.shape[0], plot_new.shape[1], 3))
for i in range(0, plot_new.shape[0]):
    for j in range(0, plot_new.shape[1]):
        data_3d[i][j] = a_color[str(int(plot_new[i][j]))][0]

plt.imshow(data_3d)

#%% calculate weights, get the heaviness per region / for all land cells / for all cells

weights = np.cos(np.deg2rad(h.lat))
weightsGPEX = np.cos(np.deg2rad(hGPEX.lat))

h_regional = h.weighted(mask_3D * weights).mean(dim=("lat", "lon")) # mean value per region
h_land     = h.weighted(land_mask * weights).mean(dim=("lat", "lon")) # 1 value for all land cells
h_all      = h.weighted(weights).mean(dim=("lat", "lon")) # 1 value for all cells

h_regional_GPEX = hGPEX.weighted(mask_3D_GPEX * weightsGPEX).mean(dim=("lat", "lon")) # mean value per region
h_land_GPEX     = hGPEX.weighted(land_mask_GPEX * weightsGPEX).mean(dim=("lat", "lon")) # 1 value for all land cells
h_all_GPEX     = hGPEX.weighted(weightsGPEX).mean(dim=("lat", "lon")) # 1 value for all cells

#%%
plot_heaviness_region = plot_new

not_included_no = [28,44,45,46,47,48,49,50,51,52,53,54,55,56,57]
for i in range(len(not_included_no)):
    plot_heaviness_region = np.where(plot_heaviness_region==not_included_no[i], np.nan, plot_heaviness_region) 

heaviness_sort = np.sort(h_regional.h10.values)

for i in range(57):
    plot_heaviness_region = np.where(plot_heaviness_region==i, h_regional.h10.values[i], plot_heaviness_region) 

plt.figure(figsize=(12,8))
plt.imshow(plot_heaviness_region,vmin=-0.5,vmax=0.5,cmap='bwr')
plt.colorbar(orientation='horizontal')
plt.ylim(601,0)
plt.title('Heaviness H10')
plt.tight_layout()
plt.show()