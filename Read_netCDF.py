# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:11:00 2021

@author: ievdv
"""
import xarray as xr
import matplotlib.pyplot as plt

#%% Import heaviness
h = xr.open_dataset('netCDF/Heaviness_H1_H5_H10.nc')

fig, ax = plt.subplots(1,3,figsize=(20,5))
h1plot = ax[0].imshow(h.h1[:601,:],vmin=-1,vmax=2)
ax[0].set_title('H1')
plt.colorbar(h1plot,ax=ax[0],orientation='horizontal')

h5plot = ax[1].imshow(h.h5[:601,:],vmin=-1,vmax=2)
ax[1].set_title('H5')
plt.colorbar(h5plot,ax=ax[1],orientation='horizontal')

h10plot = ax[2].imshow(h.h10[:601,:],vmin=-1,vmax=2)
ax[2].set_title('H10')
plt.colorbar(h10plot,ax=ax[2],orientation='horizontal')

#%% Import temperature
temp = xr.open_dataset('netCDF/Temperature.nc')

plt.figure(figsize=(9,5))
plt.imshow(temp.T_mean[:601,:])
plt.title('T mean')
plt.colorbar(orientation='horizontal')
plt.show()

fig, ax = plt.subplots(3,4,figsize=(20,10))
for i in range(12):
    plot1 = ax[i%3,i%4].imshow(temp.T_monthly[i,:601,:],vmin=230,vmax=320)
    ax[i%3,i%4].set_title(f'T mean {i+1}')
    plt.colorbar(plot1, ax=ax[i%3,i%4],orientation='horizontal')
plt.tight_layout()
plt.show()

#%% Import elevation
elev = xr.open_dataset('netCDF/Elevation.nc')

plt.figure(figsize=(9,5))
plt.imshow(elev.Elevation[:601,:])
plt.title('Elevation')
plt.colorbar(orientation='horizontal')
plt.show()

#%% Import slope aspect
slope = xr.open_dataset('netCDF/Slope.nc')

plt.figure(figsize=(9,5))
plt.imshow(slope.Aspect[:601,:])
plt.title('Aspect of the slope')
plt.colorbar(orientation='horizontal')
plt.show()

#%% Import wind direction
wind = xr.open_dataset('netCDF/Wind.nc')

plt.figure(figsize=(9,5))
plt.imshow(wind.Phi100[:601,:])
plt.title('Wind direction at 100 m')
plt.colorbar(orientation='horizontal')
plt.show()

plt.figure(figsize=(9,5))
plt.imshow(wind.Phi10[:601,:])
plt.title('Wind direction at 10 m')
plt.colorbar(orientation='horizontal')
plt.show()

#%% Import wind slope
wind = xr.open_dataset('netCDF/Windslope.nc')

plt.figure(figsize=(9,5))
plt.imshow(wind.WindSlope[:601,:])
plt.title('Wind direction and slope aspect')
plt.colorbar(orientation='horizontal')
plt.show()

#%% Import clouds
cloud = xr.open_dataset('netCDF/Clouds.nc')

fig, ax = plt.subplots(2,2,figsize=(9,5))
hc = ax[0,0].imshow(cloud.hcc[:601,:],vmin=0,vmax=1)
ax[0,0].set_title('High cloud cover')
plt.colorbar(hc,ax=ax[0,0],orientation='horizontal')

mc = ax[0,1].imshow(cloud.mcc[:601,:],vmin=0,vmax=1)
ax[0,1].set_title('Mid cloud cover')
plt.colorbar(mc,ax=ax[0,1],orientation='horizontal')

lc = ax[1,0].imshow(cloud.lcc[:601,:],vmin=0,vmax=1)
ax[1,0].set_title('Low cloud cover')
plt.colorbar(lc,ax=ax[1,0],orientation='horizontal')

tc = ax[1,1].imshow(cloud.tcc[:601,:],vmin=0,vmax=1)
ax[1,1].set_title('Total cloud cover')
plt.colorbar(tc,ax=ax[1,1],orientation='horizontal')
plt.tight_layout()
plt.show()