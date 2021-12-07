# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:51:33 2021

@author: ievdv
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import regionmask
import netCDF4 as nc
import matplotlib
import weighted_region_mean
from scipy.stats import gaussian_kde

#%% open dataset
ERA5 = xr.open_dataset('netCDF/Heaviness_H1_H5_H10.nc')
GPEX = xr.open_dataset('heaviness/GPEX.nc')
TEMP = xr.open_dataset('netCDF/Temperature.nc')
PREC = xr.open_dataset('netCDF/Precipitation.nc')
ELEV = xr.open_dataset('netCDF/Elevation.nc')
SLOPE = xr.open_dataset('netCDF/Slope.nc')
WIND = xr.open_dataset('netCDF/Wind.nc')
WINDSLOPE = xr.open_dataset('netCDF/Windslope.nc')
CLOUDS = xr.open_dataset('netCDF/Clouds.nc')

h1_era5 = ERA5.h1
h5_era5 = ERA5.h5
h10_era5 = ERA5.h10
h_gpex = GPEX.mev_heaviness[:,:,3]

#%% all IPCC reference regions 
ar6_all = regionmask.defined_regions.ar6.all
mask_3D = ar6_all.mask_3D(h1_era5)
mask_3D_gpex = ar6_all.mask_3D(h_gpex)

#%% land mask
land_110 = regionmask.defined_regions.natural_earth.land_110
land_mask = land_110.mask_3D(h1_era5)
land_mask_gpex = land_110.mask_3D(h_gpex)

#%% quick plot
ax = plt.axes(projection=ccrs.Robinson())
land_mask.sel(region=0).plot.imshow(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
plt.title('era5')
plt.show()

ax = plt.axes(projection=ccrs.Robinson())
land_mask_gpex.sel(region=0).plot.imshow(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
plt.title('gpex')
plt.show()

#%% calculate weights, get the heaviness per region / for all land cells / for all cells
h1_era5_mean = weighted_region_mean.Weighted_Region(h1_era5, mask_3D, land_mask, ar6_all)
h5_era5_mean = weighted_region_mean.Weighted_Region(h5_era5, mask_3D, land_mask, ar6_all)
h10_era5_mean = weighted_region_mean.Weighted_Region(h10_era5, mask_3D, land_mask, ar6_all)
h_gpex_mean = weighted_region_mean.Weighted_Region(h_gpex, mask_3D_gpex, land_mask_gpex, ar6_all)

weighted_region_mean.plotting(h1_era5_mean,'h1',h5_era5_mean,'h5',h10_era5_mean,'h10',h_gpex_mean,'gpex')

#%% Temperature and heaviness
temp_mean = weighted_region_mean.Weighted_Region(TEMP.T_mean, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, temp_mean, 'h1', 'temperature', 'Temperature [k]')

#%% Precipitation and heaviness
prec_mean = weighted_region_mean.Weighted_Region(PREC.Prec_mean, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, prec_mean, 'h1', 'precipitation', 'Precipitation [m]')

#%%
tp_mean['mean'] = temp_mean['mean']/(np.mean(temp_mean['mean']))**0.5 + 1/(prec_mean['mean']/(np.mean(prec_mean['mean'])))**0.5
tp_mean['abbrevs'] = temp_mean['abbrevs']
weighted_region_mean.plotting_control(h1_era5_mean, tp_mean, 'h1', r'$\sqrt{T}$ + 1/$\sqrt{P}$', 'Prec & Temp [m]')

#%% Elevation and heaviness
elev_mean = weighted_region_mean.Weighted_Region(ELEV.Elevation, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, elev_mean, 'h1', 'elevation', 'Elevation [m]')

#%% Slope and heaviness
slope_mean = weighted_region_mean.Weighted_Region(SLOPE.Aspect, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, slope_mean, 'h1', 'Slope', 'Aspect [deg]')

#%% Wind and heaviness
wind_mean = weighted_region_mean.Weighted_Region(WIND.Phi100, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, wind_mean, 'h1', 'Wind', 'Direction [deg]')

#%% Windslope and heaviness
windslope_mean = weighted_region_mean.Weighted_Region(WINDSLOPE.WindSlope, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, windslope_mean, 'h1', 'Wind and Slope', 'rho [-]')

#%% Clouds and heaviness
clouds_mean = weighted_region_mean.Weighted_Region(CLOUDS.tcc, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, clouds_mean, 'h1', 'tcc', 'Total Cloud Cover [-]')


clouds_mean_hcc = weighted_region_mean.Weighted_Region(CLOUDS.hcc, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, clouds_mean_hcc, 'h1', 'hcc', 'High Cloud Cover [-]')


clouds_mean_mcc = weighted_region_mean.Weighted_Region(CLOUDS.mcc, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, clouds_mean_mcc, 'h1', 'mcc', 'Middle Cloud Cover [-]')


clouds_mean_lcc = weighted_region_mean.Weighted_Region(CLOUDS.lcc, mask_3D, land_mask, ar6_all)

weighted_region_mean.plotting_control(h1_era5_mean, clouds_mean_lcc, 'h1', 'lcc', 'Low Cloud Cover [-]')


#%% Temperature scatterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(temp_mean, TEMP.T_mean, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Temperature [K]')

#%% Precipitation scatterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(prec_mean, PREC.Prec_mean, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Precipitation [m]')

#%% Elevation catterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(elev_mean, ELEV.Elevation, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Elevation [m]')

#%% Slope catterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(slope_mean, SLOPE.Aspect, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Slope [degrees]')

#%% Wind catterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(wind_mean, WIND.Phi100, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Wind direction [degrees]')

#%% Slopewind catterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(windslope_mean, WINDSLOPE.WindSlope, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'rho [-]')

#%% Clouds catterplots with regional data and density scatterplots of all the data
weighted_region_mean.plotting_scatter(clouds_mean, CLOUDS.tcc, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Total Cloud Cover [-]')
weighted_region_mean.plotting_scatter(clouds_mean_hcc, CLOUDS.hcc, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'High Cloud Cover [-]')
weighted_region_mean.plotting_scatter(clouds_mean_mcc, CLOUDS.mcc, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Middle Cloud Cover [-]')
weighted_region_mean.plotting_scatter(clouds_mean_lcc, CLOUDS.lcc, h1_era5, h1_era5_mean, h5_era5_mean, h10_era5_mean, 'Low Cloud Cover [-]')

#%% Temperature scatter and fit
weighted_region_mean.fit_scatter(temp_mean,h1_era5_mean,'Temperature','K')

#%% Precipitation scatter and fit
weighted_region_mean.fit_scatter(prec_mean,h1_era5_mean,'Precipitation','m')

#%% Elevation scatter and fit
weighted_region_mean.fit_scatter(elev_mean,h1_era5_mean,'Elevation','m')

#%% Slope scatter and fit
weighted_region_mean.fit_scatter(slope_mean,h1_era5_mean,'Slope','degrees')

#%% Wind scatter and fit
weighted_region_mean.fit_scatter(wind_mean,h1_era5_mean,'Wind','degrees')

#%% Slopewind scatter and fit
weighted_region_mean.fit_scatter(windslope_mean,h1_era5_mean,'rho','-')

#%% Clouds scatter and fit
weighted_region_mean.fit_scatter(clouds_mean,h1_era5_mean,'Total cloud cover','-')
weighted_region_mean.fit_scatter(clouds_mean_hcc,h1_era5_mean,'High cloud cover','-')
weighted_region_mean.fit_scatter(clouds_mean_mcc,h1_era5_mean,'Medium cloud cover','-')
weighted_region_mean.fit_scatter(clouds_mean_lcc,h1_era5_mean,'Low cloud cover','-')
