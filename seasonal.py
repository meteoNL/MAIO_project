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

#note: here S9 (old) is used, because the complete S9 includes velocities from Roderiks record, but this contains hourly records for the most recent years and therefore there are way too many observations.
names=['S4','S5','SHR','S6','S7','S8','S9 (old)','S10'] #names of sites
start=2006 #start year
end=2019 #end year plus one
num=46 #number of annual bins
ndays=365 #number of days per year
disctime=np.linspace(0,1,num) #discretize time

#empty dictionary
values={}
#read data velocity per site
for i in range(len(names)):
    name=names[i]+'.csv'
    array=np.genfromtxt(name, delimiter=';')
    array[:,0]=array[:,0]%1 #only part of the year is read (year number left out/subtracted)
    values[name[:-4]]=array[:,:4] #first four rows relevant (others contain positions and elevation)

#key indicates site    
for key in values:
    for i in range(len(values[key][:,0])):
        timeindices=np.arange(len(disctime)) #index array, used in selection below
        k=timeindices[np.abs(disctime[:]-values[key][i,0])==np.min(np.abs(disctime[:]-values[key][i,0]))] #minimalize difference between time value of velocity sample from site data and discretized time values and then select the index number at which this minimum occurs 
        values[key][i,0]=disctime[k] #replace time value of the velocity record by its now obtained discretized value

n=len(names)
#create empty array with at index 1 all discretized time values and subsequently mean for seven sites, maximum for these sites and minimum for these sites of velocity belonging to this discretized time
data=np.zeros((len(disctime),3*n+1))
data[:,0]=disctime
for key in values:
    #loop over sites/keys first (above)
    #and loop over discretized time array (below)
    for k in range(len(disctime)-1):
        
        #get proper index to put data of maximum, minimum and mean in
        for m in range(len(names)):
            if names[m]==key:
                j=m
        j+=1
        
        #now select wherever time is equal to the discretized time value and calculate max, mean, min of the velocity 
        data[k,j]=values[key][:,3][values[key][:,0]==disctime[k]].mean()
        data[k,j+n]=values[key][:,3][values[key][:,0]==disctime[k]].max()
        data[k,j+2*n]=values[key][:,3][values[key][:,0]==disctime[k]].min()
        
#visualize result
pl.figure(figsize=(12,8))
#x-axis label is day of year
xaxis=np.linspace(0,ndays,13,dtype=np.int32)
pl.title('Annual velocity cycle at all sites, mean, maximum and minimum')
#plot maximum, minimum and mean
for i in range(len(names)):
    pl.plot(data[:-1,0],data[:-1,i+1],label=names[i][:3],c=(0.14*i,1.-0.14*i,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,n+i+1],ls=':',c=(0.14*i,1.-0.14*i,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,2*n+i+1],ls='--',c=(0.14*i,1.-0.14*i,1.-0.14*i))
#improve plot
pl.grid()
pl.xticks(np.linspace(0,1,13),xaxis)
pl.xlabel('Day of year')
pl.ylabel('Velocity (m/yr)')
pl.legend()
pl.show()

#plot with mean values
pl.figure(figsize=(12,8))
#pl.title('Annual velocity cycle averaged over available years')
for i in range(len(names)):
   # pl.scatter(data[:-1,0],data[:-1,i+1],label=names[i],c=(0.14*i,0.14,1.-0.14*i))
    pl.plot(data[:-1,0],data[:-1,i+1],label=names[i][:3],c=(0.14*i,1.-0.14*i,1.-0.14*i))
pl.grid()
#add text: sitenames
day=153/365.
pl.text(day,190,'SHR')
pl.text(day,160,'S5')
pl.text(day,140,'S4')
pl.text(day,120,'S7')
pl.text(day,100,'S9')
pl.text(day,90,'S6')
pl.text(day,75,'S8')
pl.text(day,50,'S10')

#improve plot
pl.xticks(np.linspace(0,1,13),xaxis)
pl.xlabel('Day of year')
pl.ylabel('Velocity (m/yr)')
pl.ylim(0,275)
#pl.legend()
pl.show()