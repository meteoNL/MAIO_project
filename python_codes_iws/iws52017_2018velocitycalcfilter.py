# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 11:52:31 2018

@author: fam
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")
mpl.rcParams.update({'font.size': 21})
r=6.37e6
ndays=365.
navg=672#nbumber of hours to average signal 
nhours=24.
nmin=60.
ndeg=360.
radconv=ndeg/(2.*np.pi)
add_date=np.array([0,31,59,90,120,151,181,212,243,273,304,334])
crit_std=10.0

fn='iws52017_2018'
#fn='iws_S6'
fn=fn+'.txt'#andere file '.'
f=open(fn)

dataset=np.array([])


for line in f: 
    
    #split a record in segments containing date, spatial information, etc
    line=line.split('\t')
#    line=line.split(',')
    
    if line[0] !='\n':
        #conversion to date: first ordinary years
        if int(float(line[0]))%4!=0:
            date=int(float(line[0]))+add_date[int(float(line[1]))-1]/ndays+(int(float(line[2]))-1)/ndays+int(float(line[3]))/(nhours*ndays)+int(float(line[4]))/(nhours*ndays*nmin)
        #leap year case
        else:
            date=int(float(line[0]))+add_date[int(float(line[1]))-1]/(ndays+1)+(int(float(line[2]))-1)/(ndays+1)+int(int(float(line[3])))/(nhours*(ndays+1))+int(int(float(line[4])))/(nhours*(ndays+1)*nmin)
            if int(float(line[1])) > 2.5:
                date+=1./(ndays+1)
        #end of date conversion
        
        #get position
        lon=float(line[46])
        lat=float(line[47])
        z=float(line[48])
        
        #add to dataset
        dataset=np.append(dataset,[date,lon,lat,z])

f.close()

#reshape dataset to t,x,y,z-rows
dataset=dataset.reshape((len(dataset)/4,4))
time=dataset[:,0]
timecom=dataset[:,0]
lons=dataset[:,1]
lats=dataset[:,2]
height=dataset[:,3]

#height=height[height!=0.]
#lons=lons[height!=0.]
#lats=lats[height!=0.]
#timecom=timecom[height!=0.]

#mi=0
#f=0
for y in range (0,len(height)): #in deze file zijn er punten handmatig uit gehaald, timecom is de nieuwe lijst met tijden waarbij tijden van verwijderde punten eruit zijn. evt kan ik de data interpoleren naar de oude tijdslijst (met begin en eind eruit gehaald, voor zover dat voor hoogte is gedaan) 
    if height[y]==0.:
        height[y]=-10000.
        lons[y]=-10000.
        lats[y]=-10000.
        timecom[y]=-10000.

height=height[height!=-10000.]
lons=lons[lons!=-10000.]
lats=lats[lats!=-10000.]
timecom=timecom[timecom!=-10000.]

#%%
pl.figure(figsize=(12,8))
pl.title('raw data')
pl.plot(timecom, lons,'r')
#pl.plot(timecom, lats,'g')
pl.grid()
pl.ylabel('lon/lat (degree)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('raw data')
#pl.plot(timecom, lons,'r')
pl.plot(timecom, lats,'g')
pl.grid()
pl.ylabel('lon/lat (degree)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('raw data')
pl.plot(timecom, height,'b')
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()

#%%


#hold=len(height)
#hei=height[:500]
#height=height[height>100]
#istart=hold-len(height)
istart=75
iend=11
lats=lats[istart:]
lons=lons[istart:]
timecom=timecom[istart:]
height=height[istart:]
height=height[:-iend]
lats=lats[:-iend]
lons=lons[:-iend]
timecom=timecom[:-iend]
print(len(timecom), len(lats), len(lons), len(height))



def filter_dataset():
    mean_height=np.mean(height)
    std_height=np.std(height)

    for i in range(len(height)):
        height_i=height[i]
        if np.abs(height_i-mean_height)/std_height > crit_std: 
            height_i=0.5*height[(i-1)]+0.5*height[(i+1)]
            lat_i=0.5*lats[(i-1)]+0.5*lats[(i+1)]
            lon_i=0.5*lons[(i-1)]+0.5*lons[(i+1)]
            lons[i],lats[i],height[i]=lon_i,lat_i,height_i
    mean_height=np.mean(height)
    std_height=np.std(height)


filter_dataset()

meanlonssc=np.mean(lons)
meanlatssc=np.mean(lats)
meanheightsc=np.mean(height)
standheightsc=np.std(height)
heightscat1=np.subtract(height,meanheightsc)
heightscat=np.divide(heightscat1,standheightsc)
stdlinepl=np.ones([len(heightscat)])*10.
stdlinemin=np.ones([len(heightscat)])*(-10.)





pl.figure(figsize=(12,8))
pl.title('raw data')
pl.plot(timecom, lons,'r')
#pl.plot(timecom, lats,'g')
pl.grid()
pl.ylabel('lon/lat (degree)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('raw data')
#pl.plot(timecom, lons,'r')
pl.plot(timecom, lats,'g')
pl.grid()
pl.ylabel('lon/lat (degree)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('raw data')
pl.plot(timecom, height,'b')
pl.grid()
pl.ylabel('Height (m)')
pl.xlabel('Time (yr)')
pl.show()

#%%
meanheightsc2=np.mean(height)
standheightsc2=np.std(height)
heightscat12=np.subtract(height,meanheightsc2)
heightscat2=np.divide(heightscat12,standheightsc2)
pl.figure()
pl.plot(heightscat,'bo')
pl.plot(heightscat2,'ko')
pl.plot(stdlinepl,'r')
pl.plot(stdlinemin,'r')
#%%
time=timecom
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
        if dt < (1.25*navg/(2.*ndays*nhours)) and dt > 0.00:
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

dataaveraged1=np.column_stack((tavg,lonsavg,latsavg,zavg))#in text file eerst omschrijven naar strings
dataaveraged2=np.column_stack((timeval,uval,vval,velval))
#open("dataaveragedlonslats"+str(fn)+".npy", "a")
#np.save("dataaveragedlonslats"+str(fn)+".npy", dataaveraged1)
#open("dataaveragedvelocities"+str(fn)+".npy", "a")
#np.save("dataaveragedlonslats"+str(fn)+".npy", dataaveraged2)


np.save('tlonlatz'+str(fn[:-4])+'.npy', dataaveraged1) 
np.save('tuvtotu'+str(fn[:-4])+'.npy', dataaveraged2)




#text_file = open("dataaveragedlonslats"+str(fn)+".txt", "a")
#dataaveraged1.save("dataaveragedlonslats"+str(fn)+".npy")#np.load(naam)
##text_file.write(dataaveraged1)
#text_file.close()
#text_file = open("dataaveragedvelocities"+str(fn)+".txt", "a")
##text_file.write(dataaveraged2)
#text_file.close()