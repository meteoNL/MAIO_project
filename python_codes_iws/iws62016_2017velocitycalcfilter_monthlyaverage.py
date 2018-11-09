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

fn='iws62016_2017'
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
        
np.save('datasetS6test1.npy', dataset)        
       

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
lat1=lats
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
#height=height[height>100]
#istart=hold-len(height)
istart=22
iend=52
height=height[istart:]
lats=lats[istart:]
lons=lons[istart:]
timecom=timecom[istart:]
height=height[:-iend]
lats=lats[:-iend]
lons=lons[:-iend]
timecom=timecom[:-iend]
print(len(timecom), len(lats), len(lons), len(height))


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
#velocity correlations
start=2015 #starting year balance data
end=2018 #ending year balance data
min_length=int(1190) #*2 because measuring freq is per half an hour
#correlates summer and annual velocity, given the site name
#initiate storage array
ymbvel=np.zeros((12,3))
ymbvels=np.zeros((12,3))

lonsavg=np.zeros([(end-start)*12])
latsavg=np.zeros([(end-start)*12])
zavg=np.zeros([(end-start)*12])
tavg=np.zeros([(end-start)*12])
startdate=[-14.,17.,45.,76.,106.,137.,167.,198.,229.,259.,290.,320.]
startdates=(1./365.)*np.array(startdate)
enddate=[14.,45.,73.,104.,134.,165.,195.,226.,257.,287.,318.,348.]
enddates=(1./365.)*np.array(enddate)
j=0
#loop over all years
for y in range(start,end):
    #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
#    stm=0.
#    enm=0.085

    for i in range (0,len(startdates)):
        yr=y+startdates[i]
        yrend=y+enddates[i]
        largerthlo=lons[timecom>yr]
        largerthla=lats[timecom>yr]
        tlargerth=timecom[timecom>yr]
        largerthhei=height[timecom>yr]
        smallerthlo=largerthlo[tlargerth<yrend]
        smallerthla=largerthla[tlargerth<yrend]
        smallerthhei=largerthhei[tlargerth<yrend]
        tsmallerth=tlargerth[tlargerth<yrend]
        print (len(smallerthlo))
        if len(smallerthlo)>min_length:
            meanlo=np.mean(smallerthlo)
            meanla=np.mean(smallerthla)
            meantime=np.mean(tsmallerth)
            meanheight=np.mean(smallerthhei)
            #add year number from balance array and mass balance as well as mean velocity to array with results
            lonsavg[j]=meanlo
            latsavg[j]=meanla
            zavg[j]=meanheight
            tavg[j]=meantime
            j=j+1
            
lonsavg=lonsavg[tavg>0.] 
latsavg=latsavg[tavg>0.] 
zavg=zavg[tavg>0.] 
tavg=tavg[tavg>0.]            
            
            
            
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
pl.savefig("averagedvelocity"+str(navg)+".png")
f.close()

dataaveraged1=np.column_stack((tavg,lonsavg,latsavg,zavg))#in text file eerst omschrijven naar strings
dataaveraged2=np.column_stack((timeval,uval,vval,velval))
#open("dataaveragedlonslats"+str(fn)+".npy", "a")
#np.save("dataaveragedlonslats"+str(fn)+".npy", dataaveraged1)
#open("dataaveragedvelocities"+str(fn)+".npy", "a")
#np.save("dataaveragedlonslats"+str(fn)+".npy", dataaveraged2)


np.save('monthlytlonlatz'+str(fn[:-4])+'.npy', dataaveraged1) 
np.save('monthlytuvtotu'+str(fn[:-4])+'.npy', dataaveraged2)




#text_file = open("dataaveragedlonslats"+str(fn)+".txt", "a")
#dataaveraged1.save("dataaveragedlonslats"+str(fn)+".npy")#np.load(naam)
##text_file.write(dataaveraged1)
#text_file.close()
#text_file = open("dataaveragedvelocities"+str(fn)+".txt", "a")
##text_file.write(dataaveraged2)
#text_file.close()