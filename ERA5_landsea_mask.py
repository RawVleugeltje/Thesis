# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:32:12 2021

@author: ievdv
"""
import numpy as np

def ERA5grid(mask,lons,data):
    rlon = len(lons)
    rlonhalf = int(rlon/2)

    maskeast = mask[:,rlonhalf:]
    maskwest = mask[:,:rlonhalf]
    maskr    = np.hstack((maskeast, maskwest))
    
    data = data*maskr
    
    lons_new = lons[lons>=180]-360
    lons_old = lons[lons<180]
    lons = np.concatenate((lons_old,lons_new))
    
    data_new = np.zeros_like(data)
    lons_new = np.zeros_like(lons)
    
    for i in range(len(lons)):
        if lons[i] < 0:
            data_new[:,i-720] = data[:,i]
            lons_new[i-720] = lons[i]
        else:
            data_new[:,i+720] = data[:,i]
            lons_new[i+720] = lons[i]
    
    return data_new, lons_new