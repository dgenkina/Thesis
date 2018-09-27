# -*- coding: utf-8 -*-
"""
Created on Thu Nov 09 14:43:26 2017

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy


datafile=np.load('09Nov2017_files_373-402.npz')
kick = 'neg'

time=datafile['tlist']
fractions=np.zeros((time.size,3))
fractions[:,0]=datafile['fractionM']
fractions[:,1]=datafile['fraction0']
fractions[:,2]=datafile['fractionP']
qlist=datafile['qlist']
tRecoils=time*Erecoil/hbar

figure=plt.figure()
pan=figure.add_subplot(111)
pan.plot(time,qlist,'bo')
pan.set_xlabel('kickDelay')
pan.set_ylabel('measured q')
#pan.set_title('odtkick = -0.2')

figure2=plt.figure()
pan=figure2.add_subplot(111)
pan.plot(time,fractions[:,0],'ro')
pan.plot(time,fractions[:,1],'go')
pan.plot(time,fractions[:,2],'bo')
pan.set_xlabel('kickDelay')
pan.set_ylabel('measured fraction')
#pan.set_title('odtkick = -0.2')


sort=np.argsort(time)
time=time[sort]
qlist=qlist[sort]

for i in range(qlist.size-1):
    if kick=='pos':
        if qlist[i+1]>qlist[i]+0.8:
            qlist[i+1]=qlist[i+1]-2.0
    if kick=='neg':
        if qlist[i+1]<qlist[i]-1.0:
            qlist[i+1]=qlist[i+1]+2.0

        
A,B,dA,dB=lineFit(time,qlist,'kick time',r'quasimomentum [$k_L$]')