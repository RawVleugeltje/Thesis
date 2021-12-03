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
plot_new = np.zeros_like(h.h1)

for i in range(56):
    a = mask_3D.isel(region=i).values
    for q in range(601):
        for p in range(1440):
            if a[q,p] == True:
                plot_new[q,p] = 1*i
                
plot_lat = land_mask.lat.values
plot_lon = land_mask.lon.values

#%% Make netcdf file
try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass
ncfile = nc.Dataset('plot_new4.nc',mode='w',format='NETCDF4_CLASSIC') 
print(ncfile)

lat_dim = ncfile.createDimension('lat', len(plot_lat))     # latitude axis
lon_dim = ncfile.createDimension('lon', len(plot_lon))    # longitude axis

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees'
lon.long_name = 'longitude'
# Define a 3D variable to hold the data
plot = ncfile.createVariable('h1',np.float64,('lat','lon')) # note: unlimited dimension is leftmost
plot.units = '-' # degrees Kelvin
plot.standard_name = 'Regions IPCC' # this is a CF standard name

nlats = len(lat_dim); nlons = len(lon_dim)
lat[:] = plot_lat
lon[:] = plot_lon
plot[:,:] = plot_new

ncfile.close(); print('Dataset is closed!')

#%%
plot_new_nc = xr.open_dataset('plot_new4.nc')

#%%
from collections import defaultdict

a_color = defaultdict(list)
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
for i in range(56):
     a_color[str(i)].append(np.array(matplotlib.colors.to_rgb(cmaplist[4*i])))
     
data_3d = np.ndarray(shape=(plot_new_nc.h1.values.shape[0], plot_new_nc.h1.values.shape[1], 3))
for i in range(0, plot_new_nc.h1.values.shape[0]):
    for j in range(0, plot_new_nc.h1.values.shape[1]):
        data_3d[i][j] = a_color[str(int(plot_new_nc.h1.values[i][j]))][0]

plt.imshow(data_3d)
        
#%%
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

ax = plt.axes(projection=ccrs.PlateCarree())
plot_new_nc.h1.plot.imshow(ax=ax,transform=ccrs.PlateCarree(),cmap=cmap, add_colorbar=False)
plt.show()

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
plot_heaviness_region = plot_new_nc.h1.values

not_included_no = [28,46,47,50,48,51,53,54,55,49,52,45,44]
for i in range(len(not_included_no)):
    plot_heaviness_region = np.where(plot_heaviness_region==not_included_no[i], np.nan, plot_heaviness_region) 

heaviness_sort = np.sort(h_regional.h1.values)

for i in range(56):
    plot_heaviness_region = np.where(plot_heaviness_region==i, h_regional.h1.values[i], plot_heaviness_region) 

plt.imshow(plot_heaviness_region)
plt.show()

a_color = defaultdict(list)
cmap = plt.cm.viridis
cmaplist = [cmap(i) for i in range(cmap.N)]
for i in range(55):
      a_color[str(heaviness_sort[i])].append(np.array(matplotlib.colors.to_rgb(cmaplist[4*i])))
a_color[str(heaviness_sort[55])].append(np.array(matplotlib.colors.to_rgb('white')))
     
data_3d = np.ndarray(shape=(plot_heaviness_region.shape[0], plot_heaviness_region.shape[1], 3))
for i in range(0, plot_heaviness_region.shape[0]):
    for j in range(0, plot_heaviness_region.shape[1]):
        data_3d[i][j] = a_color[str(plot_heaviness_region[i][j])][0]
        
plt.figure(figsize=(10,5))
plt.imshow(data_3d)
# for i in range(100, 101): #plot_heaviness_region.shape[0]):
#     for j in range(100, 101): #plot_heaviness_region.shape[1]):
#         c = plot_heaviness_region[i,j]
#         plt.text(i, j, str(np.round(c,2)), va='center', ha='center')
#         print(c)
plt.show()
#%% create a dataframe with the mean heaviness per region

list_region_no  = []
list_abbrevs    = []
list_names      = []
list_means      = []

for rno in range(len(mask_3D.region.values)+2):
    if rno < 56:
        mean = h_regional.h1.isel(region=rno).values
    elif rno == 56:
        mean = h_land.h1.isel(region=0).values
    elif rno == 58:
        mean = h_all.values
    list_means.append(mean)
    list_region_no.append(rno)
    if rno < 56:
        list_abbrevs.append(ar6_all[rno].abbrev)
        list_names.append(ar6_all[rno].name)
    elif rno == 56:
        list_abbrevs.append('GLD')
        list_names.append('Global.land')
    elif rno == 57:
        list_abbrevs.append('GAL')
        list_names.append('Global.all')

df = pd.DataFrame()
df['mean']    = list_means
df['region_no'] = list_region_no
df['abbrevs']   = list_abbrevs
df['names']     = list_names

list_region_no  = []
list_abbrevs    = []
list_names      = []
list_means      = []

for rno in range(len(mask_3D_GPEX.region.values)+2):
    if rno < 56:
        mean = h_regional_GPEX.isel(region=rno).values
    elif rno == 56:
        mean = h_land_GPEX.isel(region=0).values
    elif rno == 58:
        mean = h_all_GPEX.values
    list_means.append(mean)
    list_region_no.append(rno)
    if rno < 56:
        list_abbrevs.append(ar6_all[rno].abbrev)
        list_names.append(ar6_all[rno].name)
    elif rno == 56:
        list_abbrevs.append('GLD')
        list_names.append('Global.land')
    elif rno == 57:
        list_abbrevs.append('GAL')
        list_names.append('Global.all')

dfGPEX = pd.DataFrame()
dfGPEX['mean']    = list_means
dfGPEX['region_no'] = list_region_no
dfGPEX['abbrevs']   = list_abbrevs
dfGPEX['names']     = list_names

#%% Exlude ocean and antartica
not_included = ['RAR','ARO','NPO','NAO','EPO','EAO','ARS','BOB','EIO','SPO','SAO','SIO','SOO','WAN','EAN']
not_included_no = [28,46,47,50,48,51,53,54,55,49,52,45,44]
included = list(set(list_abbrevs) - set(not_included))
df_new = pd.DataFrame()
a = []
b = []
c = []
for i in range(len(df)):
    for j in range(len(included)):
        if df['abbrevs'][i] == included[j]:
            a.append(df['mean'][i])
            b.append(dfGPEX['mean'][i])
            c.append(df['abbrevs'][i])

df_new['mean'] = a
df_new['mean GPEX'] = b
df_new['abbrevs'] = c