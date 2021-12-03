# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:20:47 2021

@author: ievdv
"""
import xarray as xr
import glob
import numpy as np
import LandSea_mask
import ERA5_landsea_mask

#%% Read data
list_of_paths = glob.glob('ERA5/era5_total_precipitation_*', recursive=True)

a = len(list_of_paths)
total_yearsum = []

for i in range(a):
    total_yearsum.append(xr.open_dataset(list_of_paths[i]).tp[0,:,:].values)
    
data = xr.open_dataset(list_of_paths[0])
lats = data.latitude
lons = data.longitude

#%% Precipitation land sea mask

yearly_prec = np.mean(total_yearsum, axis=0)

mask = LandSea_mask.create_mask_array(lats, lons, shapefile='continent')
yearly_prec_new, lons_new = ERA5_landsea_mask.ERA5grid(mask, lons, yearly_prec)
