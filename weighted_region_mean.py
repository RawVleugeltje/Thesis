# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 12:11:01 2021

@author: ievdv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def Weighted_Region(data, mask_3D, land_mask, ar6_all):
    weights = np.cos(np.deg2rad(data.lat))
    
    regional = data.weighted(mask_3D * weights).mean(dim=("lat", "lon")) # mean value per region
    land     = data.weighted(land_mask * weights).mean(dim=("lat", "lon")) # 1 value for all land cells
    all_cells= data.weighted(weights).mean(dim=("lat", "lon")) # 1 value for all cells

    list_region_no  = []
    list_abbrevs    = []
    list_names      = []
    list_means      = []

    for rno in range(len(mask_3D.region.values)):
        if rno < 56:
            mean = regional.isel(region=rno).values
        elif rno == 56:
            mean = land.isel(region=0).values
        elif rno == 58:
            mean = all_cells.values
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

    not_included = ['RAR','ARO','NPO','NAO','EPO','EAO','ARS','BOB','EIO','SPO','SAO','SIO','SOO','WAN','EAN']
    not_included_no = [28,46,47,50,48,51,53,54,55,49,52,45,44]
    included = list(set(list_abbrevs) - set(not_included))
    
    list_region_no  = []
    list_abbrevs    = []
    list_names      = []
    list_means      = []
    for i in range(len(df)):
        for j in range(len(included)):
            if df['abbrevs'][i] == included[j]:
                list_region_no.append(df['region_no'][i])
                list_abbrevs.append(df['abbrevs'][i])
                list_names.append(df['names'][i])
                list_means.append(df['mean'][i])
        
    df_new = pd.DataFrame()
    df_new['mean']    = list_means
    df_new['region_no'] = list_region_no
    df_new['abbrevs']   = list_abbrevs
    df_new['names']     = list_names
    
    return df_new

def plotting(data,label,data2=np.array([0]),label2=None,data3=np.array([0]),label3=None,data4=np.array([0]),label4=None):
    plt.figure(figsize=(12,3))
    plt.plot(data['abbrevs'],data['mean'],'ro',label=f'{label}')
    plt.plot(data['abbrevs'],data['mean'],'r')
    if data2.shape[0] != 1:
        plt.plot(data2['abbrevs'],data2['mean'],'bo',label=f'{label2}')
        plt.plot(data2['abbrevs'],data2['mean'],'b')
    if data3.shape[0] != 1:
        plt.plot(data3['abbrevs'],data3['mean'],'go',label=f'{label3}')
        plt.plot(data3['abbrevs'],data3['mean'],'g')
    if data4.shape[0] != 1:
        plt.plot(data4['abbrevs'],data4['mean'],'mo',label=f'{label4}')
        plt.plot(data4['abbrevs'],data4['mean'],'m')
    plt.legend()
    plt.xticks(rotation=90)
    plt.grid()
    plt.xlabel('IPCC region')
    plt.ylabel('Mean heaviness [-]')
    plt.show()
    
    return

def plotting_control(heaviness,control,label1,label2,yaxis):
    fig,ax = plt.subplots(figsize=(12,3))
    lns1 = ax.plot(heaviness['abbrevs'],heaviness['mean'],'ro',label=f'{label1}')
    ax.plot(heaviness['abbrevs'],heaviness['mean'],'r')
    ax.set_ylabel('Heaviness [-]')
    plt.xticks(rotation=90)
    ax2 = ax.twinx()
    lns2 = ax2.plot(control['abbrevs'],control['mean'],'bo',label=f'{label2}')
    ax2.plot(control['abbrevs'],control['mean'],'b')
    ax2.set_ylabel(f'{yaxis}')

    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='upper right')
    plt.show()
    
    return

def plotting_scatter(data_mean,data,h1,h1_mean,h5_mean,h10_mean,xlabel):
    plt.scatter(data_mean['mean'],h1_mean['mean'],label='h1')
    plt.scatter(data_mean['mean'],h5_mean['mean'],label='h5')
    plt.scatter(data_mean['mean'],h10_mean['mean'],label='h10')
    plt.xlabel(f'{xlabel}')
    plt.ylabel('Heaviness [-]')
    plt.grid()
    plt.legend()
    plt.show()

    plt.scatter(data,h1)
    plt.scatter(data_mean['mean'],h1_mean['mean'])
    plt.xlabel(f'{xlabel}')
    plt.ylabel('Heaviness [-]')
    plt.grid()
    plt.show()
    
    return

def fit_scatter(data_mean,h1,data_label,unit):
    X = data_mean.iloc[:,0].values.reshape(-1, 1)
    Y = h1.iloc[:,0].values.reshape(-1, 1)

    fit = LinearRegression().fit(X, Y)

    x = np.linspace(X.min(),X.max(),45)
    x = x.reshape(-1,1)
    Y_pred = fit.predict(x)
    plt.scatter(X, Y)
    plt.plot(x, Y_pred,'r--')
    plt.title(f'{data_label} --- R$^{2}$ = {fit.score(X,Y):.3f}')
    plt.xlabel(f'{data_label} [{unit}]')
    plt.ylabel('Heaviness [-]')
    plt.grid()
    plt.show()
    
    return fit.coef_, fit.intercept_