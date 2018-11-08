# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:48:24 2018

@author: fam
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
pl.close("all")


S5=np.load('tuvtotuiws_S5.npy')
S6=np.load('tuvtotuiws_S6.npy')
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
pl.plot(S5[:,0], S5[:,3],'b')

pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()

pl.figure(figsize=(12,8))
pl.title('S6')

pl.plot(s6[:,0], s6[:,3],'r')
pl.plot(s6y1617[:,0], s6y1617[:,3],'g')
pl.plot(s6y1718[:,0], s6y1718[:,3],'k')
pl.plot(S6[:,0], S6[:,3],'b')

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
pl.title('S9')

pl.plot(s10[:,0], s10[:,3],'r')


pl.grid()
pl.ylabel('velocity (m/s)')
pl.xlabel('Time (yr)')
pl.show()


test1=np.load('datasetS6test1.npy')
test2=np.load('datasetS6test2.npy')

print(test1==test2)
