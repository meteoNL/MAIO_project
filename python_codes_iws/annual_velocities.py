# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:48:24 2018

@author: fam
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import itertools
pl.close("all")


#S5=np.load('tuvtotuiws_S5.npy')
#S6=np.load('tuvtotuiws_S6.npy')
s5=np.load('tuvtotuiws5.npy')
s6=np.load('tuvtotuiws6.npy')
s9=np.load('tuvtotuiws9.npy')
s10=np.load('tuvtotuiws10.npy')
s5y1617=np.load('tuvtotuiws52016_2017.npy')
s5y1718=np.load('tuvtotuiws52017_2018.npy')
s6y1617=np.load('tuvtotuiws62016_2017.npy')
s6y1718=np.load('tuvtotuiws62017_2018.npy')
s9y1618=np.load('tuvtotuiws92016_2018.npy')

pl.figure(figsize=(12,8))
pl.title('S5')

pl.plot(s5[:,0], s5[:,3],'r')
pl.plot(s5y1617[:,0], s5y1617[:,3],'g')
pl.plot(s5y1718[:,0], s5y1718[:,3],'k')
#pl.plot(S5[:,0], S5[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S6')

pl.plot(s6[:,0], s6[:,3],'r')
pl.plot(s6y1617[:,0], s6y1617[:,3],'g')
pl.plot(s6y1718[:,0], s6y1718[:,3],'k')
#pl.plot(S6[:,0], S6[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S9')

pl.plot(s9[:,0], s9[:,3],'r')
pl.plot(s9y1618[:,0], s9y1618[:,3],'g')
  

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S10')

pl.plot(s10[:,0], s10[:,3],'r')


pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

S5vels=np.array([])
S5time=np.array([])
S5vels=np.append(S5vels,np.array([s5[:,3],s5y1617[:,3],s5y1718[:,3]]))
S5time=np.append(S5time,np.array([s5[:,0],s5y1617[:,0],s5y1718[:,0]]))


#S6vels=itertools.chain(s6[:,3], s6y1617[:,3], s6y1718[:,3])
#S6time=itertools.chain(s6[:,0], s6y1617[:,0], s6y1718[:,0])



S6vels=np.array([])
S6time=np.array([])
x2=s6[:,3]
xx=s6y1617[:,3]
z=s6y1718[:,3]
x1=s6[:,0]
xx1=s6y1617[:,0]
z1=s6y1718[:,0]
S6vels=np.concatenate([x2, xx, z])
S6time=np.concatenate([x1, xx1, z1])


S9vels=np.array([])
S9time=np.array([])

x9=s9[:,3]
xx9=s9y1618[:,3]

x99=s9[:,0]
xx99=s9y1618[:,0]

S9vels=np.concatenate([x9, xx9])
S9time=np.concatenate([x99, xx99])

#S9vels=np.append(S9vels,[s9[:,3],s9y1618[:,3]])
#S9time=np.append(S9time,[s9[:,0],s9y1618[:,0]])

S10vels=np.array([])
S10time=np.array([])
S10vels=np.append(S10vels,[s10[:,3]])
S10time=np.append(S10time,[s10[:,0]])


#%%
#velocity correlations
start=2016 #starting year balance data
end=2018 #ending year balance data
min_length=11*2+1

#correlates summer and annual velocity, given the site name
#initiate storage array
ymbvel=np.zeros((3,2))

#loop over all years
for y in range(start,end):
    #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
    yr=y-1./3.
    largerth=S5vels[S5time>yr]
    tlargerth=S5time[S5time>yr]
    smallerth=largerth[tlargerth<yr+1.]
    tsmallerth=tlargerth[tlargerth<yr+1.]
    if len(smallerth)>min_length:
        meanvel=np.mean(smallerth)
        meantime=np.mean(tsmallerth)
        
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel[y-start,0]=meantime
        ymbvel[y-start,1]=meanvel
        
    
#remove values that were not calculated and are therefore still initiated zeros        
ymbvel=ymbvel[ymbvel[:,1]>0]
    
title25='S5 annual velocity iws data'
pl.figure(figsize=(12,8))
pl.title(title25)
pl.plot(ymbvel[:,0],ymbvel[:,1])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()

print(ymbvel)
np.save('iwsannualS5.npy', ymbvel)

#%%
#velocity correlations
start=2016 #starting year balance data
end=2018 #ending year balance data
min_length=11*2+1

#correlates summer and annual velocity, given the site name
#initiate storage array
ymbvel6=np.zeros((3,2))

#loop over all years
for y in range(start,end):
    #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
    yr=y-1./3.
    largerth6=S6vels[S6time>yr]
    tlargerth6=S6time[S6time>yr]
    smallerth6=largerth6[tlargerth6<yr+1.]
    tsmallerth6=tlargerth6[tlargerth6<yr+1.]
    if len(smallerth6)>min_length:
        meanvel6=np.mean(smallerth6)
        meantime6=np.mean(tsmallerth6)
        
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel6[y-start,0]=meantime6
        ymbvel6[y-start,1]=meanvel6
        
    
#remove values that were not calculated and are therefore still initiated zeros        
ymbvel6=ymbvel6[ymbvel6[:,1]>0]
    
title26='S6 annual velocity iws data'
pl.figure(figsize=(12,8))
pl.title(title26)
pl.plot(ymbvel6[:,0],ymbvel6[:,1])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()
print(ymbvel6)

np.save('iwsannualS6.npy', ymbvel6)

#%%
#velocity correlations
start=2016 #starting year balance data
end=2018 #ending year balance data
min_length=11*2+1

#correlates summer and annual velocity, given the site name
#initiate storage array
ymbvel9=np.zeros((3,2))

#loop over all years
for y in range(start,end):
    #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
    yr=y-1./3.
    largerth9=S9vels[S9time>yr]
    tlargerth9=S9time[S9time>yr]
    smallerth9=largerth9[tlargerth9<yr+1.]
    tsmallerth9=tlargerth9[tlargerth9<yr+1.]
    if len(smallerth9)>min_length:
        meanvel9=np.mean(smallerth9)
        meantime9=np.mean(tsmallerth9)
        
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel9[y-start,0]=meantime9
        ymbvel9[y-start,1]=meanvel9
        
    
#remove values that were not calculated and are therefore still initiated zeros        
ymbvel9=ymbvel9[ymbvel9[:,1]>0]
print(ymbvel9)
    
title29='S9 annual velocity iws data'
pl.figure(figsize=(12,8))
pl.title(title29)
pl.plot(ymbvel9[:,0],ymbvel9[:,1])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()
np.save('iwsannualS9.npy', ymbvel9)

#%%
#velocity correlations
start=2016 #starting year balance data
end=2018 #ending year balance data
min_length=11*2+1

#correlates summer and annual velocity, given the site name
#initiate storage array
ymbvel10=np.zeros((3,2))

    



#loop over all years
for y in range(start,end):
    #annual velocity mean: selection of the velocity data with time smaller than or larger than 1 september value, calculate their mean and save in array for year y
    yr=y-1./3.
    largerth10=S10vels[S10time>yr]
    tlargerth10=S10time[S10time>yr]
    smallerth10=largerth10[tlargerth10<yr+1.]
    tsmallerth10=tlargerth10[tlargerth10<yr+1.]
    if len(smallerth10) !=0:
        print(len(smallerth10))
        print(len(smallerth10)>min_length)
        print(tsmallerth10)
    if len(smallerth10)>min_length:
        meanvel10=np.mean(smallerth10)
        meantime10=np.mean(tsmallerth10)
        
        #add year number from balance array and mass balance as well as mean velocity to array with results
        ymbvel10[y-start,0]=meantime10
        ymbvel10[y-start,1]=meanvel10
        
    
#remove values that were not calculated and are therefore still initiated zeros        
ymbvel10=ymbvel10[ymbvel10[:,1]>0]
    
title210='S10 annual velocity iws data'
pl.figure(figsize=(12,8))
pl.title(title210)
pl.plot(ymbvel10[:,0],ymbvel10[:,1])
pl.grid()
pl.ylabel('Velocity (m/yr)')
pl.xlabel('Time (yr)')
pl.show()
print(ymbvel10)

np.save('iwsannualS10.npy', ymbvel10)

#test1=np.load('datasetS6test1.npy')
#test2=np.load('datasetS6test2.npy')
#
#print(np.max(test1-test2))
