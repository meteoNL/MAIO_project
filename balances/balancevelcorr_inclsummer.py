# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:45:41 2018

@author: edward, inge
"""

#imports
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl

#getting ready for plotting and saving data
pl.close("all")
mpl.rcParams.update({'font.size': 21})
data=np.zeros((26,8))
bal=np.zeros(np.shape(data))

#read balancedata: time series set
data[:,0]=np.genfromtxt('meltimes.txt', delimiter='\t')[:-2]

#names of locations 
names=['S4','S5','S6','S7','S8','S9','SHR']
counter=1
min_length=85
min_length_s=23

#read balance data
for i in names:
    fn='cumbal'+i+'.txt'
    data[:,counter]=np.genfromtxt(fn, delimiter='\t')
    counter+=1
    
#add time and first values for first year to new balance array
bal[:,0]=data[:,0]
bal[0,:]=data[0,:]

#do differential to calculate subsequent mass balance for eacch year
bal[1:,1:]=data[1:,1:]-data[:-1,1:]
pl.plot(bal[:,0],bal[:,1:])

#visualization of mass balances
pl.title('Balances (mwe/year)')
pl.xlabel('Time (year)')
pl.ylabel('Balances (mwe/year)')
pl.grid(True)
#pl.savefig("balances.png")  
pl.show()
  
#velocity correlations
start=2005 #starting year balance data
end=2016 #ending year balance data
def corr_vel(site):
    #correlates summer and annual velocity, given the site name
    #initiate storage array
    ymbvel=np.zeros((12,3))
    ymbvels=np.zeros((12,3))
    
    #name file and import its data with velocities for the respective site with time as 0th and velocity as 3th column
    nm=site+'.csv'
    for i in range(len(names)):
        if names[i]==site:
            ind=(i+1)
    vel_data=np.genfromtxt(nm, delimiter=';')
    vel_time=vel_data[:,0]
    vel=vel_data[:,3]
    
    #loop over all years
    for y in range(start,end):
        #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
        yr=y-1./3
        largerth=vel[vel_time>yr]
        tlargerth=vel_time[vel_time>yr]
        smallerth=largerth[tlargerth<yr+1.]
        if len(smallerth)>min_length:
            meanvel=np.mean(smallerth)
            
            #add year number from balance array and mass balance as well as mean velocity to array with results
            ymbvel[y-start,0]=bal[(y-end),0]
            ymbvel[y-start,1]=meanvel
            ymbvel[y-start,2]=bal[(y-end),ind]
            
        #summer velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
        ysmin=y-0.59#0.41-1.00 is day 150 approx. 30 May
        
        #shouldn't this be September 1st = -1./3
        ysplus=y-0.32
        largerths=vel[vel_time>ysmin]
        tlargerths=vel_time[vel_time>ysmin]
        smallerths=largerths[tlargerths<ysplus]
        if len(smallerths)>min_length_s:
            meanvels=np.mean(smallerths)
            ymbvels[y-start,0]=bal[(y-end),0]
            print(bal[(y-end),0],yr,yr+1.)
            ymbvels[y-start,1]=meanvels
            ymbvels[y-start,2]=bal[(y-end),ind]
        
    #remove values that were not calculated and are therefore still initiated zeros        
    ymbvel=ymbvel[ymbvel[:,1]>0]
    ymbvels=ymbvels[ymbvels[:,1]>0]
    
    #visualize velocity and balance relations in scatterplot
    pl.figure()
    pl.scatter(ymbvel[:,1],ymbvel[:,2])
    pl.show() 
    pl.figure()
    pl.scatter(ymbvels[:,1],ymbvels[:,2], c='g')
    pl.show() 
    return ymbvel, ymbvels

#do the correlation, now for S7
corr_S7,corr_S7_summer=corr_vel('S7')
print(corr_S7)
S6_dat=corr_vel('S6')

#do all correlations and save in dictionary
correlations={}
for key in names:
    correlations[key]=corr_vel(key)

#correlation function for different time lags
def crosscor(lijst1,lijst2,timelist):
    
    #empty array to save data
    crosscorrlist=np.array([])
    cumcorr=0.
    
    #loop over any lags that doesn't have less than 3 observations (because correlation cann not be calculated with 2 obs)
    for lag in range (0,int(len(lijst1)-1)):
        if len(lijst2[lag:])>2:
            
            #calculate correlation coefficient
            crosscorr=np.corrcoef(lijst1[:(len(lijst1)-lag)],lijst2[lag:])[0,1]
            cumcorr=cumcorr+crosscorr
            lagyears=timelist[lag]-timelist[0]
            
            #append it to the array that saves data
            crosscorrlist=np.append(crosscorrlist,[lagyears,crosscorr,cumcorr])
    return crosscorrlist



##### Some visualzations #### ---- Vanaf hier ben ik het overzicht helemaal kwijt + ik had een vraag waarom we eigenlijk cumulatieve correlaties hadden



crosscorrlist=crosscor(corr_S7[:,2],corr_S7[:,1],corr_S7[:,0])
crosscorrlist=crosscorrlist.reshape((len(crosscorrlist)/3,3))
pl.figure()
pl.plot(crosscorrlist[:,0],crosscorrlist[:,1],'b')
pl.plot(crosscorrlist[:,0],crosscorrlist[:,2],'r')
pl.title('Cross correlation',fontsize=22)
pl.ylabel('Cross correlation',fontsize=16)
pl.xlabel('Lag',fontsize=16)
pl.grid(True)

crosscorrlists=crosscor(corr_S7_summer[:,2],corr_S7_summer[:,1],corr_S7_summer[:,0])
crosscorrlists=crosscorrlists.reshape((len(crosscorrlists)/3,3))
pl.figure()
pl.plot(crosscorrlists[:,0],crosscorrlists[:,1],'b')
pl.plot(crosscorrlists[:,0],crosscorrlists[:,2],'r')
pl.title('Cross correlation summer velocities',fontsize=22)
pl.ylabel('Cross correlation',fontsize=16)
pl.xlabel('Lag',fontsize=16)
pl.grid(True)