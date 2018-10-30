# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:00:39 2018

@author: Inge, Edward
"""

#imports
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl

#set plotting settings
pl.close("all")
mpl.rcParams.update({'font.size': 21})

names=['SHR','S4','S5','S6','S7','S8','S9 (old)','S10'] #names of sites
start=2006 #start year
end=2019 #end year plus one
num=46
ndays=365
disctime=np.linspace(0,1,num)#[:-1]

#empty dictionary
values={}
#read data velocity per site
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    array[:,0]=array[:,0]%1 #only part of the year is read (year left out)
    values[name[:-4]]=array[:,:4]
    
for key in values:
    for i in range(len(values[key][:,0])):
        timeindices=np.arange(len(disctime))
        k=timeindices[np.abs(disctime[:]-values[key][i,0])==np.min(np.abs(disctime[:]-values[key][i,0]))]
        values[key][i,0]=disctime[k]

n=len(names)
data=np.zeros((len(disctime),3*n+1))
data[:,0]=disctime
for key in values:
    for k in range(len(disctime)-1):
        for m in range(len(names)):
            if names[m]==key:
                j=m
        j+=1
        data[k,j]=values[key][:,3][values[key][:,0]==disctime[k]].mean()
        data[k,j+n]=values[key][:,3][values[key][:,0]==disctime[k]].max()
        data[k,j+2*n]=values[key][:,3][values[key][:,0]==disctime[k]].min()
pl.figure(figsize=(12,8))
xaxis=np.linspace(0,ndays,13,dtype=np.int32)
pl.title('Annual velocity cycle at all sites, mean, maximum and minimum')
for i in range(len(names)):
    pl.plot(data[:-1,0],data[:-1,i+1],label=names[i][:3],c=(0.14*i,1.-0.14*i,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,n+i+1],ls=':',c=(0.14*i,1.-0.14*i,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,2*n+i+1],ls='--',c=(0.14*i,1.-0.14*i,1.-0.14*i))
pl.grid()
pl.xticks(np.linspace(0,1,13),xaxis)
pl.xlabel('Day of year')
pl.ylabel('Velocity (m/yr)')
pl.legend()
pl.show()

pl.figure(figsize=(12,8))
pl.title('Annual velocity cycle averaged over available years')
for i in range(len(names)):
   # pl.scatter(data[:-1,0],data[:-1,i+1],label=names[i],c=(0.14*i,0.14,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,i+1],label=names[i][:3],c=(0.14*i,1.-0.14*i,1.-0.14*i))
pl.grid()
pl.xticks(np.linspace(0,1,13),xaxis)
pl.xlabel('Day of year')
pl.ylabel('Velocity (m/yr)')
pl.ylim(0,275)
pl.legend()
pl.show()