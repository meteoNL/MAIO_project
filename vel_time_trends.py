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

names=['S4','S5','SHR','S6','S7','S8','S9','S10'] #names of sites
start=2006 #start year
end=2019 #end year plus one
min_length=20.5 #length of observations per year or season criterion, in number of 96 hours intervals 
season='summer'
seasonbgn=5./12
seasonend=8./12
pl.figure(figsize=(12,8))

#empty dictionary
values={}
pl.figure(figsize=(12,8))
#read data velocity per site
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    values[name[:-4]]=array
    #switch anything below on to plot all velocity time series (which has too many scatter points)
#    pl.scatter(array[:,0],array[:,3],c=(0.14*i,1.-0.14*i,1.-0.14*i),s=2,label=names[i])    
#pl.title('Velocity series (96h average) at all sites')
#pl.xlabel('Time (yr)')
#pl.ylabel('Velocity (m/yr)')
#pl.ylim(0,400)
#pl.grid()
#pl.xlim(2002,2019)
#pl.legend(loc=2)
#pl.show()        

#linear regression function with significance from scipy and print result
def linreg(time,vel):
    lin_time_trend=st.linregress(time,vel)
    print(key, lin_time_trend)
    return key, lin_time_trend

print(season+' trends')
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
      
    #read separate IWS data from their files, depending on the season
    if key == 'S5' or key == 'S6':
        fn='iws'+season+key+'.npy'
        if season == 'annual':
            new_vel[-3:-1]=np.load(fn)[:2,1]
            new_time[-3:-1]=np.array([2016.,2017.])
        elif key == 'S6' and season == 'winter':
            new_vel[-2]=np.load(fn)[0,1]
            new_time[-2]=np.array([2017.])
        elif key == 'S6' and season == 'summer':
            new_vel[-3]=np.load(fn)[0,1]
            new_time[-3]=np.array([2016.])
    #remove zeros
    new_vel=new_vel[new_time>0]
    new_time=new_time[new_time>0]
    
    #plot result
    pl.scatter(new_time,new_vel,c=(0.14*j,1.-0.14*j,1.-0.14*j))
    pl.plot(new_time,new_vel,c=(0.14*j,1.-0.14*j,1.-0.14*j),label=key)
    
    #calculate statistical time trend significance using scipy linalg regression function
    linreg(new_time,new_vel)
    j+=1
            
#pl.title(season+' average velocity series')
#arrays with y coorinate of label positions
summerlabels=np.array([96,148,160,74,126,82,102,56])
annuallabels=np.array([83,120,112,40,104,75,96,56])
winterlabels=np.array([78,112,86,40,103,70,95,56])

#make the plot nicer
pl.xlabel('Time (yr)')
pl.ylabel('Velocity (m/yr)')
pl.ylim(0,180)
pl.grid()
pl.xlim(2006,2019)
year=2017.5
#choice of correct y-coordinates, season dependent
if season=='winter':
    labelling=winterlabels
elif season=='annual':
    labelling=annuallabels
elif season=='summer':
    labelling=summerlabels
#add text of the label
for i in range(len(names)):
    pl.text(year,labelling[i],str(names[i]))
pl.show()
