# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:48:24 2018

@author: fam
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")


#S5=np.load('tuvtotuiws_S5.npy')
#S6=np.load('tuvtotuiws_S6.npy')
s5=np.load('monthlytlonlatziws5.npy')
s6=np.load('monthlytlonlatziws6.npy')
s9=np.load('monthlytlonlatziws9.npy')
s10=np.load('monthlytlonlatziws10.npy')
s5y1617=np.load('monthlytlonlatziws52016_2017.npy')
s5y1718=np.load('monthlytlonlatziws52017_2018.npy')
s6y1617=np.load('monthlytlonlatziws62016_2017.npy')
s6y1718=np.load('monthlytlonlatziws62017_2018.npy')
s9y1618=np.load('monthlytlonlatziws92016_2018.npy')

#%%
pl.figure(figsize=(12,8))
pl.title('S5 lons')

pl.plot(s5[:,0], s5[:,1],'r')
pl.plot(s5y1617[:,0], s5y1617[:,1],'g')
pl.plot(s5y1718[:,0], s5y1718[:,1],'k')
#pl.plot(S5[:,0], S5[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S5 lats')

pl.plot(s5[:,0], s5[:,2],'r')
pl.plot(s5y1617[:,0], s5y1617[:,2],'g')
pl.plot(s5y1718[:,0], s5y1718[:,2],'k')
#pl.plot(S5[:,0], S5[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

#%%
pl.figure(figsize=(12,8))
pl.title('S6')

pl.plot(s6[:,0], s6[:,2],'r')
pl.plot(s6y1617[:,0], s6y1617[:,2],'g')
pl.plot(s6y1718[:,0], s6y1718[:,2],'k')
#pl.plot(S6[:,0], S6[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S6')

pl.plot(s6[:,0], s6[:,1],'r')
pl.plot(s6y1617[:,0], s6y1617[:,1],'g')
pl.plot(s6y1718[:,0], s6y1718[:,1],'k')
#pl.plot(S6[:,0], S6[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

#%%

pl.figure(figsize=(12,8))
pl.title('S9 lats')

pl.plot(s9[:,0], s9[:,1],'r')
pl.plot(s9y1618[:,0], s9y1618[:,1],'g')
  

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S9 lats')

pl.plot(s9[:,0], s9[:,2],'r')
pl.plot(s9y1618[:,0], s9y1618[:,2],'g')
  

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()
#%%
pl.figure(figsize=(12,8))
pl.title('S10')

pl.plot(s10[:,0], s10[:,2],'r')


pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

S5lons=np.array([])
S5lats=np.array([])
S5time=np.array([])
#S5lons=np.concatenate([s5[:,1],s5y1617[:,1],s5y1718[:,1]])
#S5lats=np.concatenate([s5[:,2],s5y1617[:,2],s5y1718[:,2]])
#S5time=np.concatenate([s5[:,0],s5y1617[:,0],s5y1718[:,0]])

avlonsvelS51=(s5y1617[-1,1]-s5y1617[-2,1])/(s5y1617[-1,0]-s5y1617[-2,0])
avlonsvelS52=(s5y1718[1,1]-s5y1718[0,1])/(s5y1718[1,0]-s5y1718[0,0])
avlonsvelS5=np.mean([avlonsvelS51,avlonsvelS52])
gapS5lons=avlonsvelS5*(s5y1718[0,0]-s5y1617[-1,0])
correctionS5lons=(s5y1718[0,1]-s5y1617[-1,1])-gapS5lons
print(gapS5lons)

for i in range(len(s5y1718[:,1])):
    s5y1718[i,1] -= correctionS5lons
  
avlatsvelS51=(s5y1617[-1,2]-s5y1617[-2,2])/(s5y1617[-1,0]-s5y1617[-2,0])
avlatsvelS52=(s5y1718[1,2]-s5y1718[0,2])/(s5y1718[1,0]-s5y1718[0,0])
avlatsvelS5=np.mean([avlatsvelS51,avlatsvelS52])
gapS5lats=avlatsvelS5*(s5y1718[0,0]-s5y1617[-1,0])
correctionS5lats=(s5y1718[0,2]-s5y1617[-1,2])-gapS5lats
print(gapS5lats)

for i in range(len(s5y1718[:,2])):
  s5y1718[i,2] = s5y1718[i,2]-correctionS5lats

S5lons=np.concatenate([s5[:,1],s5y1617[:,1],s5y1718[:,1]])
S5lats=np.concatenate([s5[:,2],s5y1617[:,2],s5y1718[:,2]])
S5time=np.concatenate([s5[:,0],s5y1617[:,0],s5y1718[:,0]])

#S6vels=itertools.chain(s6[:,3], s6y1617[:,3], s6y1718[:,3])
#S6time=itertools.chain(s6[:,0], s6y1617[:,0], s6y1718[:,0])

#%%
ndeg=360.
r=6.37e6 #radius of earth in m
correction_factor=2*np.pi*r/ndeg
ndeg=360. #number of degrees in a circle
radconv=ndeg/(2.*np.pi) #conversion between degrees and radians

def velocitycalculator (lonbegin,lonend,latbegin,latend,timebegin,timeend):

            
            #x-differential, y-differential and time differential, re-expressed in m and yr
    dx=(lonend-lonbegin)*correction_factor*np.cos(latbegin/radconv)
    dy=(latend-latbegin)*correction_factor
    dt=timeend-timebegin
    timetotvel=0.5*(timeend+timebegin)
            
    #calculate velocity components and absolute velocity in m/yr
    u=dx/dt
    v=dy/dt
    velocity=(u**2.+v**2.)**(0.5)   
    return timetotvel,velocity
#%%


S6lons=np.array([])
S6lats=np.array([])
S6time=np.array([])
#%%
print('gap before:')
print(s6y1617[0,2]-s6[-1,2])
avlonsvelS6a1=(s6[-1,1]-s6[-2,1])/(s6[-1,0]-s6[-2,0])
avlonsvelS6a2=(s6y1617[1,1]-s6y1617[0,1])/(s6y1617[1,0]-s6y1617[0,0])
avlonsvelS6a=np.mean([avlonsvelS6a1,avlonsvelS6a2])
print('mean velocity:')
print(avlonsvelS6a)
gapS6alons=avlonsvelS6a*(s6y1617[0,0]-s6[-1,0])
correctionS6alons=(s6y1617[0,1]-s6[-1,1])-gapS6alons
print(gapS6alons)

for i in range(len(s6y1617[:,1])):
    s6y1617[i,1] -= correctionS6alons
  
avlatsvelS6a1=(s6[-1,2]-s6[-2,2])/(s6[-1,0]-s6[-2,0])
avlatsvelS6a2=(s6y1617[1,2]-s6y1617[0,2])/(s6y1617[1,0]-s6y1617[0,0])
avlatsvelS6a=np.mean([avlatsvelS6a1,avlatsvelS6a2])
gapS6alats=avlatsvelS6a*(s6y1617[0,0]-s6[-1,0])
#gapS6alats=avlatsvelS6a*(2016.75-2016.66)
correctionS6alats=(s6y1617[0,2]-s6[-1,2])-gapS6alats
print(gapS6alats)

for i in range(len(s6y1617[:,2])):
  s6y1617[i,2] = s6y1617[i,2]-correctionS6alats

print('gap after:')
print(s6y1617[0,2]-s6[-1,2])

velocity_after=(s6y1617[0,1]-s6[-1,1])/(s6y1617[0,0]-s6[-1,0])
print('velocity after:')
print(velocity_after)

#----------------
avlonsvelS61=(s6y1617[-1,1]-s6y1617[-2,1])/(s6y1617[-1,0]-s6y1617[-2,0])
avlonsvelS62=(s6y1718[1,1]-s6y1718[0,1])/(s6y1718[1,0]-s6y1718[0,0])
avlonsvelS6=np.mean([avlonsvelS61,avlonsvelS62])
gapS6lons=avlonsvelS6*(s6y1718[0,0]-s6y1617[-1,0])
correctionS6lons=(s6y1718[0,1]-s6y1617[-1,1])-gapS6lons
print(gapS6lons)

for i in range(len(s6y1718[:,1])):
    s6y1718[i,1] -= correctionS6lons
  
avlatsvelS61=(s6y1617[-1,2]-s6y1617[-2,2])/(s6y1617[-1,0]-s6y1617[-2,0])
avlatsvelS62=(s6y1718[1,2]-s6y1718[0,2])/(s6y1718[1,0]-s6y1718[0,0])
avlatsvelS6=np.mean([avlatsvelS61,avlatsvelS62])
gapS6lats=avlatsvelS6*(s6y1718[0,0]-s6y1617[-1,0])
correctionS6lats=(s6y1718[0,2]-s6y1617[-1,2])-gapS6lats
print(gapS6lats)

for i in range(len(s6y1718[:,2])):
  s6y1718[i,2] = s6y1718[i,2]-correctionS6lats
#%%
print(s6y1718[0,1]-s6y1617[-1,1])
print(s6y1718[0,2]-s6y1617[-1,2])

S6lons=np.concatenate([s6[:,1], s6y1617[:,1], s6y1718[:,1]])
S6lats=np.concatenate([s6[:,2], s6y1617[:,2], s6y1718[:,2]])
S6time=np.concatenate([s6[:,0], s6y1617[:,0], s6y1718[:,0]])


S9vels=np.array([])
S9time=np.array([])
#_________________
avlonsvelS91=(s9[-1,1]-s9[-2,1])/(s9[-1,0]-s9[-2,0])
avlonsvelS92=(s9y1618[1,1]-s9y1618[0,1])/(s9y1618[1,0]-s9y1618[0,0])
avlonsvelS9=np.mean([avlonsvelS91,avlonsvelS92])
gapS9lons=avlonsvelS9*(s9y1618[0,0]-s9[-1,0])
correctionS9lons=(s9y1618[0,1]-s9[-1,1])-gapS9lons
print(gapS9lons)

for i in range(len(s9y1618[:,1])):
    s9y1618[i,1] -= correctionS9lons
  
avlatsvelS91=(s9[-1,2]-s9[-2,2])/(s9[-1,0]-s9[-2,0])
avlatsvelS92=(s9y1618[1,2]-s9y1618[0,2])/(s9y1618[1,0]-s9y1618[0,0])
avlatsvelS9=np.mean([avlatsvelS91,avlatsvelS92])
gapS9lats=avlatsvelS9*(s9y1618[0,0]-s9[-1,0])
correctionS9lats=(s9y1618[0,2]-s9[-1,2])-gapS9lats
print(gapS9lats)

for i in range(len(s9y1618[:,2])):
  s9y1618[i,2] = s9y1618[i,2]-correctionS9lats
  #__________________________________

S9lons=np.concatenate([s9[:,1], s9y1618[:,1]])
S9lats=np.concatenate([s9[:,2], s9y1618[:,2]])
S9time=np.concatenate([s9[:,0], s9y1618[:,0]])

#S9vels=np.append(S9vels,[s9[:,3],s9y1618[:,3]])
#S9time=np.append(S9time,[s9[:,0],s9y1618[:,0]])

S10vels=np.array([])
S10time=np.array([])
S10lons=s10[:,1]
S10lats=s10[:,2]
S10time=s10[:,0]
#%% S5 velocities
#S5winter=np.zeros((3,2))
#S5summer=np.zeros((3,2))
#S5winter1time,S5winter1=velocitycalculator (S5lons[0],S5lons[9],S5lats[0],S5lats[9],S5time[0],S5time[9])
#S5winter2time,S5winter2=velocitycalculator (S5lons[12],S5lons[21],S5lats[12],S5lats[21],S5time[12],S5time[21])
#S5winter3time,S5winter3=velocitycalculator (S5lons[24],S5lons[34],S5lats[24],S5lats[34],S5time[24],S5time[34])
#S5summer1time,S5summer1=velocitycalculator (S5lons[9],S5lons[12],S5lats[9],S5lats[12],S5time[9],S5time[12])
#S5summer2time,S5summer2=velocitycalculator (S5lons[21],S5lons[24],S5lats[21],S5lats[24],S5time[21],S5time[24])
#S5summer3time,S5summer3=velocitycalculator (S5lons[34],S5lons[37],S5lats[34],S5lats[37],S5time[34],S5time[37])
#print('S5velocities:')
#print(S5winter1,S5winter2,S5winter3,S5summer1,S5summer2,S5summer3)
#
#S5winter[0,0]=S5winter1time
#S5winter[0,1]=S5winter1
#S5winter[1,0]=S5winter2time
#S5winter[1,1]=S5winter2
#S5winter[2,0]=S5winter3time
#S5winter[2,1]=S5winter3
#S5summer[0,0]=S5summer1time
#S5summer[0,1]=S5summer1
#S5summer[1,0]=S5summer2time
#S5summer[1,1]=S5summer2
#S5summer[2,0]=S5summer3time
#S5summer[2,1]=S5summer3
#np.save('iwswinterS5.npy', S5winter) 
#np.save('iwssummerS5.npy', S5summer) 

#%% S6 velocities
S6winter=np.zeros((1,2))
S6summer=np.zeros((1,2))
S6winter1time,S6winter1=velocitycalculator (S6lons[11],S6lons[20],S6lats[11],S6lats[20],S6time[11],S6time[20])
#S6winter2time,S6winter2=velocitycalculator (S6lons[12],S6lons[21],S6lats[12],S6lats[21],S6time[12],S6time[21])
#S6winter3time,S6winter3=velocitycalculator (S6lons[24],S6lons[34],S6lats[24],S6lats[34],S6time[24],S6time[34])
S6summer1time,S6summer1=velocitycalculator (S6lons[8],S6lons[11],S6lats[8],S6lats[11],S6time[8],S6time[11])
#S6summer2time,S6summer2=velocitycalculator (S6lons[21],S6lons[24],S6lats[21],S6lats[24],S6time[21],S6time[24])
#S6summer3time,S6summer3=velocitycalculator (S6lons[34],S6lons[37],S6lats[34],S6lats[37],S6time[34],S6time[37])

print('S6velocities:')
print(S6winter1,S6summer1)

S6winter[0,0]=S6winter1time
S6winter[0,1]=S6winter1
#S6winter[1,0]=S6winter2time
#S6winter[1,1]=S6winter2
#S6winter[2,0]=S6winter3time
#S6winter[2,1]=S6winter3
S6summer[0,0]=S6summer1time
S6summer[0,1]=S6summer1
#S6summer[1,0]=S6summer2time
#S6summer[1,1]=S6summer2
#S6summer[2,0]=S6summer3time
#S6summer[2,1]=S6summer3
np.save('iwswinterS6.npy', S6winter) 
np.save('iwssummerS6.npy', S6summer)

##%% S9 velocities
#S9winter=np.zeros((3,2))
#S9summer=np.zeros((3,2))
#S9winter1time,S9winter1=velocitycalculator (S9lons[0],S9lons[9],S9lats[0],S9lats[9],S9time[0],S9time[9])
#S9winter2time,S9winter2=velocitycalculator (S9lons[12],S9lons[22],S9lats[12],S9lats[22],S9time[12],S9time[22])
#S9winter3time,S9winter3=velocitycalculator (S9lons[25],S9lons[34],S9lats[25],S9lats[34],S9time[25],S9time[34])
#S9summer1time,S9summer1=velocitycalculator (S9lons[9],S9lons[12],S9lats[9],S9lats[12],S9time[9],S9time[12])
#S9summer2time,S9summer2=velocitycalculator (S9lons[22],S9lons[25],S9lats[22],S9lats[25],S9time[22],S9time[25])
#S9summer3time,S9summer3=velocitycalculator (S9lons[34],S9lons[37],S9lats[34],S9lats[37],S9time[34],S9time[37])
#
#print('S9velocities:')
#print(S9winter1,S9winter2,S9winter3,S9summer1,S9summer2,S9summer3)
#S9winter[0,0]=S9winter1time
#S9winter[0,1]=S9winter1
#S9winter[1,0]=S9winter2time
#S9winter[1,1]=S9winter2
#S9winter[2,0]=S9winter3time
#S9winter[2,1]=S9winter3
#S9summer[0,0]=S9summer1time
#S9summer[0,1]=S9summer1
#S9summer[1,0]=S9summer2time
#S9summer[1,1]=S9summer2
#S9summer[2,0]=S9summer3time
#S9summer[2,1]=S9summer3
#np.save('iwswinterS9.npy', S9winter) 
#np.save('iwssummerS9.npy', S9summer)

#%% S10 velocities no velocities, not enough data points!!!!
#S10winter1time,S10winter1=velocitycalculator (S10lons[0],S10lons[9],S10lats[0],S10lats[9],S10time[0],S10time[9])
#S10winter2time,S10winter2=velocitycalculator (S10lons[12],S10lons[21],S10lats[12],S10lats[21],S10time[12],S10time[21])
#S10winter3time,S10winter3=velocitycalculator (S10lons[24],S10lons[33],S10lats[24],S10lats[33],S10time[24],S10time[33])
#S10summer1time,S10summer1=velocitycalculator (S10lons[9],S10lons[12],S10lats[9],S10lats[12],S10time[9],S10time[12])
#S10summer2time,S10summer2=velocitycalculator (S10lons[21],S10lons[24],S10lats[21],S10lats[24],S10time[21],S10time[24])
#
#print('S10velocities:')
#print(S10winter1,S10winter2,S10winter3,S10summer1,S10summer2)


#NO velocities,not enough datapoints!!!