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
min_length_s=21
min_length_w=66

#conversion mass balance year and calender year
bgnyr=-1./3
endyr=+2./3
bgns=+5./12
ends=+2./3
bgnw=bgnyr
endw=bgns

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
start=1994 #starting year balance data
end=2016 #ending year balance data

def period(y,begin,end,vel,time,min_len):
    #calculates mean velocity over a period of the year for year y plus or minus begin and end value, as selected from the time and velocity arrays delivered and with given required number of records available
    yr=y+begin
    ys=y+end
    largerth=vel[time>yr]
    tlargerth=time[time>yr]
    smallerth=largerth[tlargerth<ys]
    if len(smallerth)>min_len:
        meanvel=np.mean(smallerth)
    else:
        meanvel=np.nan
    return meanvel

def corr_vel(site):
    #correlates summer and annual velocity, given the site name
    #initiate storage array
    ymbvel,ymbvels,ymbvelw=np.zeros((22,3)),np.zeros((22,3)),np.zeros((22,3))
    
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
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel[y-start,0]=bal[(y-end),0]
        ymbvel[y-start,1]=period(y,bgnyr,endyr,vel,vel_time,min_length)
        ymbvel[y-start,2]=bal[(y-end),ind]

        #write balances also to winter and summer array
        ymbvels[y-start,0],ymbvelw[y-start,0]=bal[(y-end),0],bal[(y-end),0]
        ymbvels[y-start,2],ymbvelw[y-start,2]=bal[(y-end),ind],bal[(y-end),ind]
        
        #calculate seasonal mean and add to array
        ymbvels[y-start,1]=period(y,bgns,ends,vel,vel_time,min_length_s)
        ymbvelw[y-start,1]=period(y,bgnw,endw,vel,vel_time,min_length_w)
        
    #visualize velocity and balance relations in scatterplot
    pl.figure()
    pl.scatter(ymbvel[:,1],ymbvel[:,2])
    pl.show() 
    pl.figure()
    pl.scatter(ymbvels[:,1],ymbvels[:,2], c='g')
    pl.show() 
    return ymbvel, ymbvels,ymbvelw

#do all correlations and save in dictionary
correlations={}
for key in names:
    correlations[key]=corr_vel(key)

#correlation function for different time lags
def crosscor(array):
    lijst1,lijst2,timelist=array[:,2],array[:,1],array[:,0]
    crosscorrlist=np.array([])

    
    #loop over any lags that doesn't have less than 3 observations (because correlation cann not be calculated with 2 obs)
    for lag in range (0,8):
        if len(lijst2[lag:])>2:
            
            #calculate correlation coefficient
            #apply first proper selection: select x-series and y-series (mass balance and velocity) for this lag and remove nan's
            x=lijst1[:(len(lijst1)-lag)]
            y=lijst2[lag:]
            
            #remove nan's in y
            x=x[y<np.inf]
            y=y[y<np.inf]
            
            #calculate correlation coefficient
            crosscorr=np.corrcoef(x,y)[0,1]
            lagyears=timelist[lag]-timelist[0]
            
            #append correlation coefficient to the array that saves data
            crosscorrlist=np.append(crosscorrlist,[lagyears,crosscorr])
            
    #[0::2] gives lagyears, [1::2] gives crosscorr
    return crosscorrlist[0::2],crosscorrlist[1::2]

def plotting(lag,y):
    #define function for lag correlations plot
    pl.scatter(lag,y)
    pl.plot(lag,y,label=key)

#create empty dictionaries for data    
yearlagcorr={}
summerlagcorr={}
winterlagcorr={}

for key in names:
    #get data named above and calculate the correlation as function of time 
    yearcrosscorsite,summercrosscorsite,wintercrosscorsite=correlations[key][0],correlations[key][1],correlations[key][2]
    yearlagcorr[key]=crosscor(yearcrosscorsite)
    summerlagcorr[key]=crosscor(summercrosscorsite)
    winterlagcorr[key]=crosscor(wintercrosscorsite)

#plot this time series
pl.figure(figsize=(12,8))
for key in names:     
    season=', annual',plotting(yearlagcorr[key][0],yearlagcorr[key][1])
    #season=' summer',plotting(summerlagcorr[key][0],summerlagcorr[key][1])
    #season=' winter', plotting(winterlagcorr[key][0],winterlagcorr[key][1])
pl.title('Balance and velocity relation'+str(season[0]))
pl.ylabel('Correlation coefficient (-)')
pl.xlim(0,9)
pl.xlabel('Lag (yr)')
pl.grid()
pl.legend(loc=1)
pl.show()
