# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 16:05:10 2018

@author: Inge, Edward
"""
#imoports, constants
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")
mpl.rcParams.update({'font.size': 21})
r=6.37e6
ndays=365.
navg=48 #nbumber of hours to average signal 
nhours=24.
ndeg=360.
radconv=ndeg/(2.*np.pi)
add_date=np.array([0,31,59,90,120,151,181,212,243,273,304,334])
crit_std=5.0

#filename
fn='S4'
fn=fn+'.txt'
f=open(fn)

dataset=np.array([])

for line in f: 
    
    #split a record in segments containing date, spatial information, etc
    line=line.split(',')
    
    if line[0] !='\n':
        #conversion to date: first ordinary years
        if int(line[0][0:2])%4!=0:
            date=2e3+int(line[0][0:2])+add_date[int(line[0][-4:-2])-1]/ndays+(int(line[0][-2:])-1)/ndays+int(line[1][:2])/(nhours*ndays)
        #leap year case
        else:
            date=2e3+int(line[0][0:2])+add_date[int(line[0][-4:-2])-1]/(ndays+1)+(int(line[0][-2:])-1)/(ndays+1)+int(line[1][:2])/(nhours*(ndays+1))
            if int(line[0][-4:-2]) > 2.5:
                date+=1./(ndays+1)
        #end of date conversion
        
        #get position
        lon=float(line[2])
        lat=float(line[3])
        z=float(line[4])
        
        #add to dataset
        dataset=np.append(dataset,[date,lon,lat,z])

f.close()
    
#reshape dataset to t,x,y,z-rows
dataset=dataset.reshape((len(dataset)/4,4))
time=dataset[:,0]
lons=dataset[:,1]
lats=dataset[:,2]
height=dataset[:,3]

def filter_dataset():
    mean_height=np.mean(height)
    std_height=np.std(height)
    counter=0
    while np.max(abs(height-np.mean(height))/np.std(height)) > crit_std:
        for i in range(len(dataset[:,3])):
            height_i=dataset[i,3]
            if np.abs(height_i-mean_height)/std_height > crit_std:
                height_i=0.5*dataset[(i-1),3]+0.5*dataset[(i+1),3]
                lat_i=0.5*dataset[(i-1),2]+0.5*dataset[(i+1),2]
                lon_i=0.5*dataset[(i-1),1]+0.5*dataset[(i+1),1]
                dataset[i,1],dataset[i,2],dataset[i,3]=lon_i,lat_i,height_i
        mean_height=np.mean(height)
        std_height=np.std(height)
        counter+=1
    print(counter)

filter_dataset()

def compute_velocities(nstep,lons,lats,time,height):
    raw_velocities=np.array([])
    
    #calculations for each interval
    for i in range(len(lons)-nstep):
        
        #correction factor from degrees to metres (approx 111 km)
        correction_factor=2*np.pi*r/ndeg
        
        #x-differential, y-differential and time differential, re-expressed in m and yr
        dx=(lons[i+nstep]-lons[i])*correction_factor*np.cos(lats[i]/radconv)
        dy=(lats[i+nstep]-lats[i])*correction_factor
        dt=time[i+nstep]-time[i]
        time_t=0.5*(time[i+nstep]+time[i])
        
        #calculate velocity components and absolute velocity in m/yr
        u=dx/dt
        v=dy/dt
        velocity=(u**2.+v**2.)**(0.5)
        
        #put result in an array (inculding time)
        if dt < (1.25*navg/(ndays*nhours)) and dt > 0.00:
            raw_velocities=np.append(raw_velocities,np.array([time_t,u,v,velocity]))
    return raw_velocities[0::4],raw_velocities[1::4],raw_velocities[2::4],raw_velocities[3::4]

#some visualization
title1=fn[:2]+' vertical coordinate evolution in time'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(time,height)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()

#%%
def averagingcomputation():
    
    #create arrays
    lonsavg=np.array([])
    latsavg=np.array([])
    zavg=np.array([])
    tavg=np.array([])
    #dt=np.array([])
    
    for i in range(0,int(len(lons)/navg)):
        startind=i*navg
        endind=(i+1)*navg
        lonsavg=np.append(lonsavg,np.mean(lons[startind:endind]))
        latsavg=np.append(latsavg,np.mean(lats[startind:endind]))
        zavg=np.append(zavg,np.mean(height[startind:endind]))
        tavg=np.append(tavg,np.mean(time[startind:endind]))

    return lonsavg,latsavg,zavg,tavg
        
lonsavg,latsavg,zavg,tavg=averagingcomputation()
timeval,uval,vval,velval=compute_velocities(1,lonsavg,latsavg,tavg,zavg)

#some visualization
title1=fn[:2]+' vertical coordinate evolution in time, '+str(navg)+' hours average'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(tavg,zavg)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()
pl.figure(figsize=(12,8))
title2=fn[:2]+' absolute velocity evolution in time, '+str(navg)+' hours average'
pl.title(title2)
pl.plot(timeval,velval)
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()
f.close()