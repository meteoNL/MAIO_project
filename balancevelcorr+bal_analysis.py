# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:45:41 2018

@author: edward, inge
"""

#imports
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import scipy as sp
import scipy.stats as st

#getting ready for plotting and saving data
pl.close("all")
mpl.rcParams.update({'font.size': 21})
data=np.zeros((26,8))
bal=np.zeros(np.shape(data))

#read balancedata: time series set
data[:,0]=np.genfromtxt('meltimes.txt', delimiter='\t')[:-2]

#names of locations 
names=['S4','S5','SHR','S6','S7','S8','S9']
counter=1
#note that previously calculated S9 velocities by Roderik automatically pass these criteria, but they were assessed and found complete.  Define criteria
min_length=80.5
min_length_s=20.5
min_length_w=60.5

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
pl.figure(figsize=(12,8))
bal[1:,1:]=data[1:,1:]-data[:-1,1:]
for i in range(len(names)):
    pl.plot(bal[:,0],bal[:,i+1],c=(0.14*i,1.-0.14*i,1.-0.14*i),label=names[i])
    pl.scatter(bal[:,0],bal[:,i+1],c=(0.14*i,1.-0.14*i,1.-0.14*i))

#visualization of mass balances
#pl.title('Time series of balance rate')
pl.xlabel('Time (year)')
pl.ylabel('Balance rate (mwe/yr)')
pl.grid(True)
pl.legend(loc=1)
pl.xlim(1990,2025)  
pl.show()
  
#velocity correlations
start=1990 #starting year balance data
end=2019 #ending velocity balance data
endb=2016 #end balance data

def period(y,begin,end,vel,time,min_len):
    #calculates mean velocity over a period of the year for year y plus or minus begin and end value, as selected from the time and velocity arrays delivered and with given required number of records available
    yr=y+begin
    ys=y+end
    largerth=vel[time>yr]
    tlargerth=time[time>yr]
    smallerth=largerth[tlargerth<ys]
    
    #print if on the edge of criterion (if necessary)
  #  if len(smallerth) < 1.05*min_len and len(smallerth)!=0:
   #     print(len(smallerth),key,y)
    if len(smallerth)>min_len:
        meanvel=np.mean(smallerth)
    else:
        meanvel=np.nan
    return meanvel

def read_veldata(site):
    #name file and import its data with velocities for the respective site with time as 0th and velocity as 3th column
    nm=site+'.csv'
    for i in range(len(names)):
        if names[i]==site:
            ind=(i+1)
    vel_data=np.genfromtxt(nm, delimiter=';')
    return vel_data,ind


def corr_vel(site):
    #correlates summer and annual velocity, given the site name
    #initiate storage array
    ymbvel,ymbvels,ymbvelw=np.zeros((29,3)),np.zeros((29,3)),np.zeros((29,3))
    
    #read the velocity data and select the appropriate columns
    vel_data,ind=read_veldata(site)
    vel_time=vel_data[:,0]
    vel=vel_data[:,3]
    
    #loop over all years
    for y in range(start,endb):
        #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y            
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel[y-start,0]=bal[(y-endb),0]
        ymbvel[y-start,1]=period(y,bgnyr,endyr,vel,vel_time,min_length)
        ymbvel[y-start,2]=bal[(y-endb),ind]

        #write balances also to winter and summer array
        ymbvels[y-start,0],ymbvelw[y-start,0]=bal[(y-endb),0],bal[(y-endb),0]
        ymbvels[y-start,2],ymbvelw[y-start,2]=bal[(y-endb),ind],bal[(y-endb),ind]
        
        #calculate seasonal mean and add to array
        ymbvels[y-start,1]=period(y,bgns,ends,vel,vel_time,min_length_s)
        ymbvelw[y-start,1]=period(y,bgnw,endw,vel,vel_time,min_length_w)
        
    #add data from 2016 and 2017 and 2018, no balance data
    for y in range(endb,end):
        ymbvel[y-start,0],ymbvels[y-start,0],ymbvelw[y-start,0]=y+0.5,y+0.5,y+0.5
        ymbvel[y-start,1],ymbvels[y-start,1],ymbvelw[y-start,1]=period(y,bgnyr,endyr,vel,vel_time,min_length),period(y,bgns,ends,vel,vel_time,min_length_s),period(y,bgnw,endw,vel,vel_time,min_length_w)
        ymbvel[y-start,2],ymbvels[y-start,2],ymbvelw[y-start,2]=np.nan,np.nan,np.nan
    return ymbvel, ymbvels,ymbvelw

#do all correlations and save in dictionary
#keys are the sitenames
correlations={}
for key in names:
    correlations[key]=corr_vel(key)
    
##### ADD IWS DATA FROM THEIR FILEs ####
correlations['S5'][0][-3:-1,1]=np.load('iwsannualS5.npy')[:,1]
correlations['S6'][0][-3:-1,1]=np.load('iwsannualS6.npy')[:,1]
correlations['S6'][1][-3,1]=np.load('iwssummerS6.npy')[:,1]
correlations['S6'][2][-2,1]=np.load('iwswinterS6.npy')[:,1]

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
            
            #remove nan's in y and x
            x=x[y<np.inf]
            y=y[y<np.inf]
            y=y[x<np.inf]
            x=x[x<np.inf]
            
            #calculate correlation coefficient
            crosscorr=st.linregress(x,y)[2]
            lagyears=timelist[lag]-timelist[0]
            if (lag == 1 or lag == 2) and key in names[-4:]: 
                print(key,lag,st.linregress(x,y))
            
            #append correlation coefficient to the array that saves data
            crosscorrlist=np.append(crosscorrlist,[lagyears,crosscorr])
            
    #[0::2] gives lagyears, [1::2] gives crosscorr
    return crosscorrlist[0::2],crosscorrlist[1::2]

def plotting(lag,y,j):
    #define function for lag correlations plot
    pl.scatter(lag,y,c=(0.14*j,1.-0.14*j,1.-0.14*j))
    pl.plot(lag,y,label=key,c=(0.14*j,1.-0.14*j,1.-0.14*j))

def betterplot():
    #improve plotting
#    pl.title('Balance rate and velocity relation'+str(season[0]))
    pl.ylabel('Correlation coefficient (-)')
    pl.xlim(0,9)
    pl.ylim(-1,1)
    pl.xlabel('Lag (yr)')
    pl.grid()
    pl.legend(loc=1)
    pl.show()

#create empty dictionaries for data    
yearlagcorr={}
summerlagcorr={}
winterlagcorr={}

for key in names:
    #get year, mass balance and velocity averages per season/annum for each site (=key) 
    yearcrosscorsite,summercrosscorsite,wintercrosscorsite=correlations[key][0],correlations[key][1],correlations[key][2]
    #and subsequently calculate the correlation as function of time lag for the seasons
    yearlagcorr[key]=crosscor(yearcrosscorsite)
    summerlagcorr[key]=crosscor(summercrosscorsite)
    winterlagcorr[key]=crosscor(wintercrosscorsite)

#plot this time series: lag correlations found above
pl.figure(figsize=(12,8))
j=0
for key in names:     
    season=', annual',plotting(yearlagcorr[key][0],yearlagcorr[key][1],j)
    j+=1
betterplot()
pl.figure(figsize=(12,8))
j=0
for key in names:
    season=' summer',plotting(summerlagcorr[key][0],summerlagcorr[key][1],j)
    j+=1
betterplot()
pl.figure(figsize=(12,8))
j=0
for key in names:
    season=' winter', plotting(winterlagcorr[key][0],winterlagcorr[key][1],j)
    j+=1
betterplot()

#linear regression function that prints most important results and returns the values
def linreg(x,y):
    lin_time_trend=st.linregress(x,y)
    print(names[i], lin_time_trend)
    return lin_time_trend

#vertical and time mass balance relation plot: start with an array with 0's
verttrend=np.zeros(((len(bal[0,:])-1),2))
pl.figure(figsize=(12,8))
for i in range(len(bal[0,:])-1):
    site=names[i]
    
    #time trend for appropriate site
    result=linreg(bal[:,0],bal[:,i+1])
    
    #read height data (which were in the velocity results files)
    vel_data,ind=read_veldata(site)
    
    #calculate site mean height and balance and then put in array
    hgt=np.mean(vel_data[:,6][vel_data[:,6]<np.inf])
    balrate=np.mean(bal[:,i+1])
    verttrend[i,0],verttrend[i,1]=hgt,balrate
    
#plot the result: vertical relation of mass balance
for i in range(len(bal[0,:])-1):
    pl.scatter(verttrend[i,0],verttrend[i,1],s=50,label=names[i],marker='o',c=(0.14*i,1.-0.14*i,1.-0.14*i))
print('balance/height ')
#regression result and plot regression result subsequently
regresult=linreg(verttrend[:,0],verttrend[:,1])
hgtnew=np.linspace(verttrend[0,0],verttrend[-1,0],20)
pl.plot(hgtnew,regresult[0]*hgtnew+regresult[1],c='black',label='Fit')
pl.xlim(0,2000)
pl.ylim(-6,2)
pl.legend(loc=2)
pl.xlabel('Height above SL (m)')
pl.ylabel('Balance rate (mwe/yr)')
#pl.title('Balance rate gradient relation observed in K-transect')
pl.grid()
pl.show()