# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:45:41 2018

@author: fam
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")
mpl.rcParams.update({'font.size': 21})
data=np.zeros((26,8))
bal=np.zeros(np.shape(data))
data[:,0]=np.genfromtxt('meltimes.txt', delimiter='\t')[:-2]
names=['S4','S5','S6','S7','S8','S9','SHR']
counter=1
for i in names:
    fn='cumbal'+i+'.txt'
    data[:,counter]=np.genfromtxt(fn, delimiter='\t')
    counter+=1
bal[:,0]=data[:,0]
bal[0,:]=data[0,:]
bal[1:,1:]=data[1:,1:]-data[:-1,1:]
pl.plot(bal[:,0],bal[:,1:])
pl.title('Balances (mwe/year)')
pl.xlabel('Time (year)')
pl.ylabel('Balances (mwe/year)')
pl.grid(True)
#pl.savefig("balances.png")  
pl.show()
  
#velocity correlations
start=2005
end=2018
def corr_vel(site):
    ymbvel=np.zeros((14,3))
    nm=site+'.csv'
    for i in range(len(names)):
        if names[i]==site:
            ind=(i+1)
    vel_data=np.genfromtxt(nm, delimiter=';')
    vel_time=vel_data[:,0]
    vel=vel_data[:,3]
    for y in range(start,end):
        largerth=vel[vel_time>y]
        tlargerth=vel_time[vel_time>y]
        smallerth=largerth[tlargerth<y+1]
        if len(smallerth)>90:
            meanvel=np.mean(smallerth)
            ymbvel[y-start,0]=bal[y-start,0]
            ymbvel[y-start,1]=meanvel
            ymbvel[y-start,2]=bal[y-start,ind]
            
    ymbvel=ymbvel[ymbvel[:,1]>0]
    #pl.figure()
    #pl.scatter(ymbvel[:,1],ymbvel[:,2])
    #pl.show()
    print(np.corrcoef(ymbvel[:,1],ymbvel[:,2]))    
    return ymbvel
corr_vel('S7')