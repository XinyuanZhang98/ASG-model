# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:48:02 2023

@author: xinyuan zhang
"""

import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, zeros
from math import radians, cos, sin, asin, sqrt, pi
from sklearn.cluster import DBSCAN
import os
from os import listdir
from datetime import datetime
import gc
import math
from shapely.geometry import Point, Polygon, shape
import shapely.affinity

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    m = 6367000 * c
    return m 

def SDE(lon_list,lat_list):
    mean_lon=np.mean(lon_list)
    mean_lat=np.mean(lat_list)
    A=0
    B=0
    C=0
    D=0
    SDEX=0
    SDEY=0
    for i in range(len(lon_list)):
        delta_lon=lon_list[i]-mean_lon
        delta_lat=lat_list[i]-mean_lat
        #delta_lon=haversine(lon_list[i],mean_lat,mean_lon,mean_lat)
        #delta_lat=haversine(mean_lon,mean_lat,mean_lon,lat_list[i])
        SDEX+=delta_lon*delta_lon/len(lon_list)
        SDEY+=delta_lat*delta_lat/len(lat_list)
        A+=delta_lon*delta_lon-delta_lat*delta_lat
        B+=delta_lon*delta_lat
        C+=delta_lon*delta_lon
        D+=delta_lat*delta_lat
    
    alpha=np.arctan((A+sqrt(A*A+4*B*B))/(2*B))
    if math.isnan(alpha):
        alpha=pi/2
    if len(lon_list)==2:
        x=delta_lon/(2*cos(alpha))
        y=x
        return mean_lon,mean_lat,alpha,x,y
    else:
        
        x=sqrt(2*(cos(alpha)*cos(alpha)*C-2*sin(alpha)*cos(alpha)*B+sin(alpha)*sin(alpha)*D)/(len(lon_list)))
        y=sqrt(2*(sin(alpha)*sin(alpha)*C+2*sin(alpha)*cos(alpha)*B+cos(alpha)*cos(alpha)*D)/(len(lon_list)))
        
        return mean_lon,mean_lat,alpha,x,y

def get_index(location,region):
    for i in range(len(region)):
        if guangzhou_cb.iloc[i,7].contains(location):
            break
    return i

#calculate the weighted number of activity spaces
guangzhou_cb=gpd.read_file(r"J:\results\zxy\guangzhou_2km.shp")    
path=r'J:\results\zxy\\cluster_feature_guangzhou'
filenames=os.listdir(path)
n=0
space_covered=[0]*len(guangzhou_cb)
guangzhou_cb['space_cove']=space_covered
for filename in filenames:
    n+=1
    print(n)
    spath1=r'J:\results\zxy\\cluster_feature_guangzhou\\'+filename
    data=pd.read_csv(spath1, sep=',', header=0)
    cluster_list=list(data['cluster'])
    spath2=r'J:\results\zxy\\cluster_guangzhou\\'+filename
    #lon lat leave_hour
    data=pd.read_csv(spath2, sep=',', header=0)
    for c in cluster_list:
        sequence=data[data.cluster==c]
        try:
            lon,lat,alpha,a,b=SDE(list(sequence['lon']),list(sequence['lat']))
        except ValueError as e:
            continue
        circ=Point(lon,lat).buffer(1)
        ell=shapely.affinity.scale(circ,a,b,origin='centroid')
        ellr1=shapely.affinity.rotate(ell,180*alpha/pi,origin='centroid')
        if ellr1.area==0:
            continue
        index_list=guangzhou_cb[guangzhou_cb.intersects(ellr1)].index
        for m in index_list:  
            guangzhou_cb['space_cove'][m]+=guangzhou_cb['geometry'][m].intersection(ellr1).area/ellr1.area
    
guangzhou_cb.to_file(r"J:\results\zxy\\guangzhou_space.shp",driver='ESRI Shapefile', encoding='utf-8')    

