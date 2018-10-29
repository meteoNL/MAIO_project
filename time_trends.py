# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 09:29:28 2018

@author: Inge, Edward
"""

#HIERONDER WORDT DE OLD S9 GEBRUKT OMDAT DE SEASONAL CYCLE WEGING ANDERS NOG AANGEPAST MOET)

#imports
import numpy as np
import scipy as sp
import scipy.stats as st
import matplotlib.pyplot as pl
import matplotlib as mpl

#set plotting settings
pl.close("all")
mpl.rcParams.update({'font.size': 21})

names=['SHR','S4','S5','S6','S7','S8','S9','S10'] #names of sites
start=2006 #start year
end=2019 #end year plus one
min_length=81.5 #length of observations per year criterion, in number of 96 hours intervals 
season='annual'
seasonbgn=-1./3
seasonend=8./12
pl.figure(figsize=(12,8))
if season != 'annual':
    ender=-1
else: 
    ender=0

#empty dictionary
values={}
pl.figure(figsize=(12,8))
#read data velocity per site
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    values[name[:-4]]=array

    pl.scatter(array[:,0],array[:,3],c=(0.14*i,1.-0.14*i,1.-0.14*i),s=2,label=names[i])
    #pl.plot(array[:,0],array[:,3],c=(0.14*i,1.-0.14*i,1.-0.14*i),label=key)
      
pl.title('Velocity series (96h average) at all sites')
pl.xlabel('Time (yr)')
pl.ylabel('Velocity (m/yr)')
pl.ylim(0,400)
pl.grid()
pl.xlim(2002,2019)
pl.legend(loc=2)
pl.show()        

#linear regression function with significance from scipy
def linreg(time,vel):
    lin_time_trend=st.linregress(time,vel)
    print(key, lin_time_trend)
    return key, lin_time_trend
    
#96 h time series time trend
for key in values:
    time=values[key][:,0]
    vel=values[key][:,3]
    lin_time_trend=linreg(time,vel)
print(season+'trends')
pl.figure(figsize=(12,8))
#annual velocities time trend  
j=0  
for key in names:
    
    #read time and velocity array per site
    time=values[key][:,0]
    vel=values[key][:,3]
    
    #create new array for yearly avgs
    new_time=np.zeros((end-start))
    new_vel=np.zeros((end-start))
    i=0 #counter
    for y in range(start,end):
        
        #selection of velocity data per glaciological/mass balance year (??) year
        first=y+seasonbgn
        last=y+seasonend
        vel_yr=vel[time>first]
        time_yr=time[time>first]
        
        #if criterion fulfilled: calculate mean velocity and save in array with year
        if len(time_yr[time_yr<last]) > min_length:
            vel_yr=np.mean(vel_yr[time_yr<last])
            new_time[i]=y
            new_vel[i]=vel_yr
        i+=1
      #  print(key,y,vel_yr)
    if key == 'S5' or key == 'S6':
        fn='iws'+season+key+'.npy'
        if season == 'annual':
            new_vel[-3:-1]=np.load(fn)[:2,1]
            new_time[-3:-1]=np.array([2016.,2017.])
        else:
            new_vel[-3:]=np.load(fn)[:,1]
            new_time[-3:]=np.array([2016.,2017.,2018.])
    #remove zeros
    new_vel=new_vel[new_time>0]
    new_time=new_time[new_time>0]
    
    #plot result
    pl.scatter(new_time,new_vel,c=(0.14*j,1.-0.14*j,1.-0.14*j))
    pl.plot(new_time,new_vel,c=(0.14*j,1.-0.14*j,1.-0.14*j),label=key)
    
    #calculate statistical time trend significance using scipy linalg regression function
    linreg(new_time,new_vel)
    j+=1
            
pl.title(season+' average velocity series')
pl.xlabel('Time (yr)')
pl.ylabel('Velocity (m/yr)')
pl.ylim(0,180)
pl.grid()
pl.xlim(2002,2019)
pl.legend(loc=2)
pl.show()        
    
