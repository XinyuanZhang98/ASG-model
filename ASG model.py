# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 12:13:08 2023

@author: xinyuan zhang
"""
import pandas as pd
import numpy as np
import statsmodels.formula.api as sm
import matplotlib.pyplot as plt
import geopandas as gpd
from math import radians, cos, sin, asin, sqrt, pi
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


def CPC(realvalue, predictvalue):
    a=0
    b=0
    c=0
    mea=0
    rmse=0
    for i in range(len(realvalue)):
        a+=2*min(realvalue[i],predictvalue[i])
        b+=realvalue[i]+predictvalue[i]
        c+=2*min(realvalue[i],predictvalue[i])/(realvalue[i]+predictvalue[i])
        mea+=abs(realvalue[i]-predictvalue[i])/len(realvalue)
        rmse+=(realvalue[i]-predictvalue[i])**2/len(realvalue)
    return a/b , c/len(realvalue),mea,sqrt(rmse)



distance=np.zeros((len(guangzhou_pop),len(guangzhou_pop)))
for i in range(len(guangzhou_pop)):
    print(i)
    for j in range(len(guangzhou_pop)):
        distance[i,j]=haversine(guangzhou_pop['geometry'][i].centroid.x,guangzhou_pop['geometry'][i].centroid.y,
                          guangzhou_pop['geometry'][j].centroid.x,guangzhou_pop['geometry'][j].centroid.y)

'''
lon1=112+57/60
lon2=114+3/60
lat1=22+26/60
lat2=23+56/60
accuracy=2000
deltaLon=accuracy*360/(2*math.pi*6371004*math.cos((lat1+lat2)*math.pi/360))
deltaLat=accuracy*360/(2*math.pi*6371004)
'''

od=np.load(r'J:\results\zxy\\od_g.npy')
guangzhou_cb=gpd.read_file(r"J:\results\zxy\\guangzhou_space.shp")    
wa=np.zeros((len(guangzhou_cb),len(guangzhou_cb)))
wp=np.zeros((len(guangzhou_cb),len(guangzhou_cb)))
Ti=[0]*4788
for i in range(len(guangzhou_cb)):
    print(block,i)
    for j in range(len(guangzhou_cb)):
        if i!=j and od[i,j]!=0:
            if guangzhou_cb['pop'][i]!=0 and guangzhou_cb['pop'][j]!=0 and od[i,j]!=0 and guangzhou_cb['space_cove'][i]!=0 and guangzhou_cb['space_cove'][j]!=0:
                    d=d_g[i,j]
                    ox=i//84
                    oy=i%84-1
                    dx=j//84
                    dy=j%84-1
                    index_list=[]
                    if ox>=dx:
                        x_list=range(dx,min(56,ox+ox-dx))
                    else:
                        x_list=range(max(0,ox-dx+ox),dx)
                    if oy>=dy:
                        y_list=range(dy,min(83,oy+oy-dy))
                    else:
                        y_list=range(max(0,oy-dy+oy),dy)   
                    for x in x_list:
                        for y in y_list:
                            index_list+=[x*84+y]
                    for k in index_list:
                        if d_g[i,k]<=d:
                            wa[i,j]+=guangzhou_cb['space_cove'][k]
                            wp[i,j]+=guangzhou_cb['pop'][k]
            
T=[]        
i_list=[]
j_list=[]
P1=[]
P2=[]
A1=[]
A2=[]
T_P=[]
T_A=[]
PWO1P=[]
PWO1A=[]
for i in range(len(guangzhou_cb)):
    print(i)
    for j in range(len(guangzhou_cb)):
        if guangzhou_cb['pop'][i]!=0 and guangzhou_cb['pop'][j]!=0 and od[i,j]!=0 and i!=j:
            if guangzhou_cb['space_cove'][i]>0.1**10 and guangzhou_cb['space_cove'][j]>0.1**10:
                T+=[od[i,j]]
                i_list+=[i]
                j_list+=[j]
                P1+=[guangzhou_cb['pop'][i]]
                P2+=[guangzhou_cb['pop'][j]]
                A1+=[guangzhou_cb['space_cove'][i]]
                A2+=[guangzhou_cb['space_cove'][j]]
                T_P+=[od[i,j]/(guangzhou_cb['pop'][i]*guangzhou_cb['pop'][j])]
                T_A+=[od[i,j]/(guangzhou_cb['space_cove'][i]*guangzhou_cb['space_cove'][j])]            
                PWO1P+=[1/wp[i,j]]
                PWO1A+=[1/wa[i,j]]

                
data=pd.DataFrame()
data['T']=T
data['i']=i_list
data['j']=j_list
data['P1']=P1
data['P2']=P2
data['A1']=A1
data['A2']=A2
data['T_P']=T_P
data['T_A']=T_A
data['PWO1P']=PWO1P
data['PWO1A']=PWO1A
data = data.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
data=data.reset_index()
data.to_csv(r'J:\\results\zxy\\gravity_regression_guangzhou.txt',index=False)

def log_binning(datapoint,bins):
    interval = list(np.logspace(0,np.log10(len(datapoint)+1+10),bins,base = 10))
    x_list=[1]
    y_list=[datapoint[0,0]]
    for i in range(1,len(interval)-1):
        sequence=datapoint[np.where(datapoint[:,1]<=(interval[i+1]+interval[i])/2)]
        sequence=sequence[np.where(sequence[:,1]>(interval[i]+interval[i-1])/2)]
        if len(sequence)!=0:
            x_list+=[interval[i]]
            y_list+=[np.mean(sequence[:,0])]
    return x_list,y_list

data=pd.read_csv(r'J:\\results\zxy\\gravity_regression_guangzhou.txt',header=0)

gravity9=sm.ols(formula='T_P~PWO1P -1',data=data)
res9=gravity9.fit()
print(res9.summary())

gravity11=sm.ols(formula='T_A~PWO1A -1',data=data)
res11=gravity11.fit()
print(res11.summary())


T9=[]
T11=[]
for i in range(len(data)):
    if i%1000==0:
        print(i)
    T9+=[(res9.params[0]*data['PWO1P'][i])*data['P1'][i]*data['P2'][i]]
    T11+=[(res11.params[0]*data['PWO1A'][i])*data['A1'][i]*data['A2'][i]]


print(CPC(data['T'],T9),CPC(data['T'],T11))