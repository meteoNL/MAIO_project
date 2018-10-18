# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:45:41 2018

@author: edward, inge
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
end=2016
def corr_vel(site,lag):
    ymbvel=np.zeros((12,3))
    nm=site+'.csv'
    for i in range(len(names)):
        if names[i]==site:
            ind=(i+1)
    vel_data=np.genfromtxt(nm, delimiter=';')
    vel_time=vel_data[:,0]
    vel=vel_data[:,3]
    for y in range(start,end):
        yr=y-1./3
        largerth=vel[vel_time>yr]
        tlargerth=vel_time[vel_time>yr]
        smallerth=largerth[tlargerth<yr+1.]
        if len(smallerth)>85:
            meanvel=np.mean(smallerth)
            ymbvel[y-start,0]=bal[(y-end),0]
            #print(bal[(y-end),0],yr,yr+1.)
            ymbvel[y-start,1]=meanvel
            ymbvel[y-start,2]=bal[(y-end),ind]
            
    ymbvel=ymbvel[ymbvel[:,1]>0]
    pl.figure()
    pl.scatter(ymbvel[:,1],ymbvel[:,2])
    pl.show() 
    return ymbvel
corr_S7=corr_vel('S7',0)
print(corr_S7)

def crosscor(lijst1,lijst2,timelist):
#    laglist=np.zeros([int(len(lijst1)-1)])
#    crosscorrlist=np.zeros([int(len(lijst1)-1)])
#    cumcorrlist=np.zeros([int(len(lijst1)-1)])
#    cumcorrlistabs=np.zeros([int(len(lijst1)-1)])
    crosscorrlist=np.array([])
    cumcorr=0.
    for lag in range (0,int(len(lijst1)-1)):
        if len(lijst2[lag:])>2:
            crosscorr=np.corrcoef(lijst1[:(len(lijst1)-lag)],lijst2[lag:])[0,1]
            cumcorr=cumcorr+crosscorr
            lagyears=timelist[lag]-timelist[0]
            crosscorrlist=np.append(crosscorrlist,[lagyears,crosscorr,cumcorr])
#            cumcorrlist[lag]=cumcorr
#            cumcorrlistabs[lag]=cumcorrabs
#            laglist[lag]=timelist[lag]-timelist[0]
        #laglist[lag]=lag#lags kloppen niet in deze regel, omdat er waarschijnlijk grotere gaten tussen zitten
    return crosscorrlist

crosscorrlist=crosscor(corr_S7[:,2],corr_S7[:,1],corr_S7[:,0])
crosscorrlist=crosscorrlist.reshape((len(crosscorrlist)/3,3))
pl.figure()
pl.plot(crosscorrlist[:,0],crosscorrlist[:,1],'b')
pl.plot(crosscorrlist[:,0],crosscorrlist[:,2],'r')
pl.title('Cross correlation',fontsize=22)
pl.ylabel('Cross correlation',fontsize=16)
pl.xlabel('Lag',fontsize=16)
pl.grid(True)