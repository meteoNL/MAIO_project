# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 09:29:28 2018

@author: Inge, Edward
"""

#imports
import numpy as np
import scipy as sp
import scipy.stats as st
import matplotlib.pyplot as pl
import matplotlib as mpl

#set plotting settings
pl.close("all")
mpl.rcParams.update({'font.size': 21})

names=['S4','S5','S6','S7','S8','S9','S10','SHR'] #names of sites
start=2006 #start year
end=2019 #end year plus one
min_length=85 #length of observations per year criterion, in number of 96 hours intervals 

#empty dictionary
values={}
#read data velocity per site
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    values[name[:-4]]=array

#linear regression function with significance from scipy
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
    
    #read time and velocity array per site
    time=values[key][:,0]
    vel=values[key][:,3]
    
    #create new array for yearly avgs
    new_time=np.zeros((end-start))
    new_vel=np.zeros((end-start))
    i=0 #counter
    for y in range(start,end):
        
        #selection of velocity data per glaciological/mass balance year (??) year
        first=y-1./3
        last=y+2./3
        vel_yr=vel[time>first]
        time_yr=time[time>first]
        
        #if criterion fulfilled: calculate mean velocity and save in array with year
        if len(time_yr[time_yr<last]) > min_length:
            vel_yr=np.mean(vel_yr[time_yr<last])
            new_time[i]=y
            new_vel[i]=vel_yr
        i+=1
      #  print(key,y,vel_yr)
      
    #remove zeros
    new_vel=new_vel[new_time>0]
    new_time=new_time[new_time>0]
    
    #plot result
    pl.figure(figsize=(12,8))
    pl.title(key+' velocity per year')
    pl.xlabel('Time (yr)')
    pl.ylabel('Velocity (m/yr)')
    pl.plot(new_time,new_vel)
    
    #calculate statistical time trend significance using scipy linalg regression function
    linreg(new_time,new_vel)
            
        
    
