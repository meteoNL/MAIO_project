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
melttimes=np.genfromtxt('meltimes.txt', delimiter='\t')[:-2]
cumbalS4=np.genfromtxt('cumbalS4.txt', delimiter='\t')
cumbalS5=np.genfromtxt('cumbalS5.txt', delimiter='\t')
cumbalS6=np.genfromtxt('cumbalS6.txt', delimiter='\t')
cumbalS7=np.genfromtxt('cumbalS7.txt', delimiter='\t')
cumbalS8=np.genfromtxt('cumbalS8.txt', delimiter='\t')
cumbalS9=np.genfromtxt('cumbalS9.txt', delimiter='\t')
cumbalSHR=np.genfromtxt('cumbalSHR.txt', delimiter='\t')
balS4=np.empty([len(cumbalS4)])
balS5=np.empty([len(cumbalS5)])
balS6=np.empty([len(cumbalS6)])
balS7=np.empty([len(cumbalS7)])
balS8=np.empty([len(cumbalS8)])
balS9=np.empty([len(cumbalS9)])
balSHR=np.empty([len(cumbalSHR)])
balS4[0]=cumbalS4[0]
balS5[0]=cumbalS5[0]
balS6[0]=cumbalS6[0]
balS7[0]=cumbalS7[0]
balS8[0]=cumbalS8[0]
balS9[0]=cumbalS9[0]
balSHR[0]=cumbalSHR[0]

pl.figure()
pl.plot (melttimes,cumbalS4)
pl.plot (melttimes,cumbalS5)
pl.plot (melttimes,cumbalS6)
pl.plot (melttimes,cumbalS7)
pl.plot (melttimes,cumbalS8)
pl.plot (melttimes,cumbalS9)
pl.plot (melttimes,cumbalSHR)
pl.title('Cumulative balances (mwe/year)',fontsize=20)
pl.xlabel('Time (year',fontsize=14)
pl.ylabel('Cumulative balances (mwe/year)',fontsize=14)
pl.grid(True)
#pl.savefig("cumbalances.png")  



for i in range (1,len(cumbalS4)):
    balS4[i]=cumbalS4[i]-cumbalS4[i-1]
    balS5[i]=cumbalS5[i]-cumbalS5[i-1]
    balS6[i]=cumbalS6[i]-cumbalS6[i-1]
    balS7[i]=cumbalS7[i]-cumbalS7[i-1]
    balS8[i]=cumbalS8[i]-cumbalS8[i-1]
    balS9[i]=cumbalS9[i]-cumbalS9[i-1]
    balSHR[i]=cumbalSHR[i]-cumbalSHR[i-1]
    
pl.figure()
pl.plot (melttimes,balS4)
pl.plot (melttimes,balS5)
pl.plot (melttimes,balS6)
pl.plot (melttimes,balS7)
pl.plot (melttimes,balS8)
pl.plot (melttimes,balS9)
pl.plot (melttimes,balSHR)
pl.title('Balances (mwe/year)',fontsize=20)
pl.xlabel('Time (year',fontsize=14)
pl.ylabel('Balances (mwe/year)',fontsize=14)
pl.grid(True)
#pl.savefig("balances.png")  
pl.show()    
ymbvel=np.zeros((14,3))
S7=np.genfromtxt('S7.csv', delimiter=';')
S7time=S7[:,0]
S7totvel=S7[:,3]
for y in range(2005,2018):
    largerth=S7totvel[S7time>y]
    tlargerth=S7time[S7time>y]
    smallerth=largerth[tlargerth<y+1]
    if len(smallerth)>90:
        meanvel=np.mean(smallerth)
        ymbvel[y-2005,0]=melttimes[y-2005]
        ymbvel[y-2005,1]=meanvel
        ymbvel[y-2005,2]=balS7[y-2005]
        
ymbvel=ymbvel[ymbvel[:,1]>0]
pl.figure()
pl.scatter(ymbvel[:,1],ymbvel[:,2])
pl.show()
print(np.corrcoef(ymbvel[:,1],ymbvel[:,2]))