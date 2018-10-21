# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 09:29:28 2018

@author: Inge, Edward
"""

import numpy as np
import scipy as sp
import scipy.stats as st
import matplotlib.pyplot as pl
import matplotlib as mpl
pl.close("all")
mpl.rcParams.update({'font.size': 21})
names=['S4','S5','S6','S7','S8','S9','S10','SHR']
start=2006
end=2019
min_length=85

values={}
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    values[name[:-4]]=array

def linreg(time,vel):
    lin_time_trend=st.linregress(time,vel)
    print(key, lin_time_trend)
    return key, lin_time_trend
    
#96 h time series time trend
for key in values:
    time=values[key][:,0]
    vel=values[key][:,2]
    lin_time_trend=linreg(time,vel)
print('annual trends')
#annual velocities time trend    
for key in values:
    time=values[key][:,0]
    vel=values[key][:,3]
    new_time=np.zeros((end-start))
    new_vel=np.zeros((end-start))
    i=0
    for y in range(start,end):
        first=y-1./3
        last=y+2./3
        vel_yr=vel[time>first]
        time_yr=time[time>first]
        if len(time_yr[time_yr<last]) > min_length:
            vel_yr=np.mean(vel_yr[time_yr<last])
            new_time[i]=y
            new_vel[i]=vel_yr
        i+=1
      #  print(key,y,vel_yr)
    new_vel=new_vel[new_time>0]
    new_time=new_time[new_time>0]
    pl.figure(figsize=(12,8))
    pl.title(key+' velocity per year')
    pl.xlabel('Time (yr)')
    pl.ylabel('Velocity (m/yr)')
    pl.plot(new_time,new_vel)
    linreg(new_time,new_vel)
            
        
    
