# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 16:05:10 2018

@author: Edward
"""
#imoports, constants
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")
mpl.rcParams.update({'font.size': 21})
nsteps=48 #number of time steps to use for differential calculations (velocity = dx/dt, nsteps is dt in hours if data are complete)
r=6.37e6
ndays=365.
navg=48 #nbumber of hours to average signal 
nhours=24.
ndeg=360.
radconv=ndeg/(2.*np.pi)
add_date=np.array([0,31,59,90,120,151,181,212,243,273,304,334])

#filename
fn='S720082015int'
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
        raw_velocities[i,0]=0.5*(time[i+nstep]+time[i])
        raw_velocities[i,1]=u
        raw_velocities[i,2]=v
        raw_velocities[i,3]=velocity
        
    return raw_velocities

#execute velocity calculation
velocity_data=compute_velocities(nsteps,lons,lats,time,height)

#some visualization
title1=fn[:-4]+' vertical coordinate evolution in time'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(time,height)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()
pl.figure(figsize=(12,8))
title2=fn[:-4]+' absolute velocity evolution in time'
pl.title(title2)
pl.plot(velocity_data[:,0],velocity_data[:,3])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()


#%%

ampslat=np.fft.fft(lats)
indnum=int(len(lats)/navg)
ampslat[indnum:-indnum]=0.0
latsrec=np.fft.ifft(ampslat)
ampslon=np.fft.fft(lons)
ampslon[indnum:-indnum]=0.0
lonsrec=np.fft.ifft(ampslon)

velocity_data_fft=compute_velocities(1,lonsrec,latsrec,time,height)

#some visualization
title1=fn[:-4]+' vertical coordinate evolution in time'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(time,height)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()
pl.figure(figsize=(12,8))
title2=fn[:-4]+' absolute velocity evolution in time'
pl.title(title2)
pl.plot(velocity_data_fft[indnum:-indnum,0],velocity_data_fft[indnum:-indnum,3])
pl.grid()
pl.ylabel('Velocity (m/yr)')
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
velocity_data_daily_avg=compute_velocities(1,lonsavg,latsavg,tavg,zavg)

#some visualization
title1=fn[:-4]+' vertical coordinate evolution in time, '+str(navg)+' hours average'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(tavg[1:],zavg)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()
pl.figure(figsize=(12,8))
title2=fn[:-4]+' absolute velocity evolution in time, '+str(navg)+' hours average'
pl.title(title2)
pl.plot(velocity_data_daily_avg[:,0],velocity_data_daily_avg[:,3])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()

#%%
pl.figure(figsize=(12,8))
pl.plot(tavg[:-1],dt)
pl.show()