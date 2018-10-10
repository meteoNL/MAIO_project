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
#nsteps=48 #number of time steps to use for differential calculations (velocity = dx/dt, nsteps is dt in hours if data are complete)
r=6.37e6
ndays=365.
navg=48 #nbumber of hours to average signal 
nhours=24.
ndeg=360.
radconv=ndeg/(2.*np.pi)
add_date=np.array([0,31,59,90,120,151,181,212,243,273,304,334])
crit_std=5.0

#filename
fn='S420062018int'
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
    
#reshape dataset to t,x,y,z-rows
dataset=dataset.reshape((len(dataset)/4,4))
time=dataset[:,0]
lons=dataset[:,1]
lats=dataset[:,2]
height=dataset[:,3]

pl.figure()
pl.plot(time[:-2],-ndays*nhours*(time[:-2]-time[1:-1]))
pl.show()

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
    raw_velocities=np.zeros(((len(lons)-nstep),4))
    
    #calculations for each interval
    for i in range(len(lons)-nstep):
        
        #correction factor from degrees to metres (approx 111 km)
        correction_factor=2*np.pi*r/ndeg
        
        #x-differential, y-differential and time differential, re-expressed in m and yr
        dx=(lons[i+nstep]-lons[i])*correction_factor*np.cos(lats[i]/radconv)
        dy=(lats[i+nstep]-lats[i])*correction_factor
        dt=time[i+nstep]-time[i]
        
        #calculate velocity components and absolute velocity in m/yr
        u=dx/dt
        v=dy/dt
        velocity=(u**2.+v**2.)**(0.5)
        
        #put result in an array (inculding time)
        if dt < (2.*navg/(ndays*nhours)):
            raw_velocities[i,0]=0.5*(time[i+nstep]+time[i])
            raw_velocities[i,1]=u
            raw_velocities[i,2]=v
            raw_velocities[i,3]=velocity
        
    return raw_velocities[:,0],raw_velocities[:,1],raw_velocities[:,2],raw_velocities[:,3]

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
    tavg=np.array([time[0]])
    dt=np.array([])
    
    for i in range(0,int(len(lons)/navg)):
        tavgold=tavg
        startind=i*navg
        endind=(i+1)*navg
        lonsavg=np.append(lonsavg,np.mean(lons[startind:endind]))
        latsavg=np.append(latsavg,np.mean(lats[startind:endind]))
        zavg=np.append(zavg,np.mean(height[startind:endind]))
        tavg=np.append(tavg,np.mean(time[startind:endind]))
        dtnew=tavg[-1]-tavg[-2]
        dt=np.append(dt,np.array([dtnew]))
    return tavgold,lonsavg,latsavg,zavg,tavg,dt
        
tavgold,lonsavg,latsavg,zavg,tavg,dt=averagingcomputation()
timeval,uval,vval,velval=compute_velocities(1,lonsavg,latsavg,tavg,zavg)

#some visualization
title1=fn[:2]+' vertical coordinate evolution in time, '+str(navg)+' hours average'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(tavg[1:],zavg)
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

#%%
pl.figure(figsize=(12,8))
pl.plot(tavg[:-1],dt)
pl.show()

#%%
pl.figure()
pl.plot(timeval,velval)
pl.xlim(2000,2020)
pl.ylim(0,500)