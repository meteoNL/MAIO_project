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
r=6.37e6 #radius of earth in m
ndays=365. #number of days in a normal year
navg=96 #nbumber of hours to average signal 
nhours=24. #number of hours per day
ndeg=360. #number of degrees in a circle
radconv=ndeg/(2.*np.pi) #conversion between degrees and radians
add_date=np.array([0,31,59,90,120,151,181,212,243,273,304,334]) #day of year for each starting day of a month
crit_std=10.0 #critical number of standard deviations for filtering applications
crit_dt=1.25 #critical relative (to navg) time interval to accept a velocity sample

#%% read data
#filename
fn='S6'
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

#%%selection
#remove zeros for latitude and longitude
time=time[lons<0]
height=height[lons<0]
lats=lats[lons<0]
lons=lons[lons<0]
    
def filter_dataset():
    #calculate mean and std dev for filtering dataset based on height + set counter
    mean_height=np.mean(height)
    std_height=np.std(height)
    counter=0
    
    #repeated check for the need to be filtered
    while np.max(np.abs(height-mean_height)/std_height) > crit_std:
        print(counter)
        
        #pass along all samples
        for i in range(len(height)):
            height_i=height[i]
            
            #replace all position values by linear interpolation between nearest neighbours
            if np.abs(height_i-mean_height)/std_height > crit_std: 
                height_i=0.5*height[(i-1)]+0.5*height[(i+1)]
                lat_i=0.5*lats[(i-1)]+0.5*lats[(i+1)]
                lon_i=0.5*lons[(i-1)]+0.5*lons[(i+1)]
                
                #update with obtained values in array
                lons[i],lats[i],height[i]=lon_i,lat_i,height_i
        
        #update statistics subsequently
        mean_height=np.mean(height)
        std_height=np.std(height)
        counter+=1

#apply filter function on dataset
filter_dataset()

#%% velocity computation function
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
        
        #put result in an array (inculding time) if it fulfills the time step criterion
        if dt < (crit_dt*navg/(ndays*nhours)) and dt > 0.00:
            raw_velocities=np.append(raw_velocities,np.array([time_t,u,v,velocity]))
    return raw_velocities[0::4],raw_velocities[1::4],raw_velocities[2::4],raw_velocities[3::4]

#some visualization
title1=fn[:-4]+' vertical coordinate evolution in time'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(time,height)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()

#%%
def averagingcomputation():
    
    #create arrays to save values
    lonsavg=np.array([])
    latsavg=np.array([])
    zavg=np.array([])
    tavg=np.array([])
    
    for i in range(0,int(len(lons)/navg)):
        
        #get the subarrays
        startind=i*navg
        endind=(i+1)*navg
        
        #fill the new arrays with averages based on subarray
        lonsavg=np.append(lonsavg,np.mean(lons[startind:endind]))
        latsavg=np.append(latsavg,np.mean(lats[startind:endind]))
        zavg=np.append(zavg,np.mean(height[startind:endind]))
        tavg=np.append(tavg,np.mean(time[startind:endind]))

    return lonsavg,latsavg,zavg,tavg

#function below doesnt take iterative standard deviation and mean into account yet 
def postproc():
    #generate dx and dy array and their means and std to check statistical behavior
    gradlatsavg=latsavg[1:]-latsavg[:-1]
    gradlonsavg=lonsavg[1:]-lonsavg[:-1]
    meangradlat=np.mean(gradlatsavg[:])
    meangradlon=np.mean(gradlonsavg[:])
    stdgradlat=np.std(gradlatsavg[:])
    stdgradlon=np.std(gradlonsavg[:])
    
    #translatetion array + bookkeeping 
    trans_x=0.
    trans_y=0.
    counter=0
    checkarray=np.zeros(1)
    
    while len(checkarray) > 0: 
        
        print(np.max(np.abs(gradlatsavg[:]-meangradlat))/stdgradlat)
        print(np.max(np.abs(gradlonsavg[:]-meangradlon))/stdgradlon)
        indices=np.arange(len(gradlatsavg))
        checkarray=indices[(np.abs(gradlatsavg-meangradlat))/stdgradlat > crit_std]
        checkarray=np.append((indices[(np.abs(gradlonsavg-meangradlon))/stdgradlon > crit_std]),checkarray)
        checkarray+=1
        print(checkarray)
        for i in checkarray:
            #interpolate values of velocity based on nearest neighbours and use this to calculate a translation of the sensor for both lat and lon (only one neighbour in i = 0 or i = arraylength case)
            if i > 1 and i < len(lonsavg):
                trans_xcor=0.5*(lonsavg[i-2]-lonsavg[i-1]+lonsavg[i+1]-lonsavg[i])
                trans_ycor=0.5*(latsavg[i-2]-latsavg[i-1]+latsavg[i+1]-latsavg[i])
            elif i ==0:
                trans_xcor=lonsavg[i+1]-lonsavg[i]
                trans_ycor=latsavg[i+1]-latsavg[i]
            elif i == len(lonsavg):
                trans_xcor=lonsavg[i-2]-lonsavg[i-1]
                trans_ycor=latsavg[i-2]-latsavg[i-1]
            
            #do the translation on the full array after the jump based on two coordinate difference and its expected value based on interpolated velocity
            trans_x=lonsavg[i]-lonsavg[i-1]-trans_xcor
            trans_y=latsavg[i]-latsavg[i-1]-trans_ycor
            
            #translation
            lonsavg[i:]=lonsavg[i:]-trans_x
            latsavg[i:]=latsavg[i:]-trans_y
        
        #update old values for new statistics based on newly obtained arrays and do counter
        gradlatsavg=latsavg[1:]-latsavg[:-1]
        gradlonsavg=lonsavg[1:]-lonsavg[:-1]
        meangradlat=np.mean(gradlatsavg[:])
        meangradlon=np.mean(gradlonsavg[:])
        stdgradlat=np.std(gradlatsavg[:])
        stdgradlon=np.std(gradlonsavg[:])
        counter+=1
        print(counter)
    return lonsavg,latsavg

#the execution of averaging and post-processing        
lonsavg,latsavg,zavg,tavg=averagingcomputation()
lonsavg,latsavg=postproc()
timeval,uval,vval,velval=compute_velocities(1,lonsavg,latsavg,tavg,zavg)

#some visualization
title1=fn[:-4]+' vertical coordinate evolution in time, '+str(navg)+' hours average'
pl.figure(figsize=(12,8))
pl.title(title1)
pl.plot(tavg,zavg)
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()
pl.figure(figsize=(12,8))
title2=fn[:-4]+' absolute velocity evolution in time, '+str(navg)+' hours average'
pl.title(title2)
pl.plot(timeval,velval)
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()

#writing result to file (if uncommented)
sep=';'
#g=open(fn[:-4]+'.csv','a')
#for i in range(len(timeval)-1):
 #   line=str(timeval[i])+sep+str(uval[i])+sep+str(vval[i])+sep+str(velval[i])+sep+str(latsavg[i])+sep+str(lonsavg[i])+sep+str(zavg[i])+sep+str(latsavg[i+1])+sep+str(lonsavg[i+1])+sep+str(zavg[i+1])+sep+'\n'
  #  g.write(line)
#g.close()