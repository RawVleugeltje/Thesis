# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 10:22:17 2021

@author: ievdv
"""

import os
import numpy as np
from osgeo import gdal, ogr
import xarray as xr
import matplotlib.pyplot as plt

os.environ['PROJ_LIB'] = "C:\\Users\\ievdv\\anaconda3\\Library\\share\\proj"
os.environ['GDAL_DATA'] = 'C:\\Users\\ievdv\\anaconda3\\Library\\share'

def create_mask_array(lats, lons, shapefile=None, mask_value=1, ysign=-1):
    """Result: mask variable for your dataset, 0 = sea, 1 = land """
    nlats = len(lats)
    nlons = len(lons)
    
    if shapefile == None:
        mask = np.ones((nlats,nlons))
        no_land = sum(sum(mask == 1))

    else:
        ymax  = lats[-1]
        ymin  = lats[0]
        xmax  = lons[-1]
        xmin  = lons[0]
        if xmax > 181:
            lons = lons - 180
            xmin = lons[0]
            xmax = lons[-1]
        mask_value = mask_value
        
        xres = 360/nlons
        yres = 180/nlats
        geotransform=(xmin,xres,0,ysign*ymax,0,-yres)
                
        ''' aanpassen naar eigen folder ''' 
        shape_dir  = r'LandSea/continent_clip' 
        
        try:
            src_ds = ogr.Open(shape_dir + '/' + shapefile)
            src_lyr=src_ds.GetLayer()
        except:
            src_ds = ogr.Open(shape_dir + '/' + shapefile + '.shp')
            src_lyr=src_ds.GetLayer()
        
        dst_ds = gdal.GetDriverByName('MEM').Create('', nlons, nlats, 1 ,gdal.GDT_Byte)
        dst_rb = dst_ds.GetRasterBand(1)
        dst_rb.Fill(0) #initialise raster with zeros
        dst_rb.SetNoDataValue(0)
        dst_ds.SetGeoTransform(geotransform)
        
        err = gdal.RasterizeLayer(dst_ds, [mask_value], src_lyr)
        dst_ds.FlushCache()
        
        # all continent
        mask_arr=dst_ds.GetRasterBand(1).ReadAsArray()
        mask_arr[mask_arr>0] = mask_value
        if shapefile == 'continent':
            mask_arr[-1,:] = 1
        no_land    = sum(sum(mask_arr == 1))
        mask            = mask_arr.copy()
        mask            = mask.astype('float')
        mask[mask == 0] = np.nan

    return mask