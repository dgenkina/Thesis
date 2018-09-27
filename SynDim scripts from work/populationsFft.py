# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 17:00:48 2017

@author: dng5
"""


import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy


c=1064.0/790.0#0.0#
n=7
k=0.0
epsilon=0.02
U=4.4
S=2
m0=0
delta=0.097
omega=0.45

filename = '25Jul2017_files_184-213'
dataFile=np.load(filename+'.npz')
imbal=dataFile['imbalArray']
signalGood=dataFile['signalGood']
cutoff=0.5
fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
#fieldGoodArray=np.ones(dataFile['tlist'].size,dtype=bool)


toffs=0.0
tList=dataFile['tlist'][fieldGoodArray]+toffs
sort=np.argsort(tList)
tList=tList[sort]

fractionP=dataFile['fractionP'][fieldGoodArray][sort]
fraction0=dataFile['fraction0'][fieldGoodArray][sort]
fractionM=dataFile['fractionM'][fieldGoodArray][sort]


maxInd=np.int(fractionP.size)
fractions=np.append(np.append(fractionM[:maxInd],fraction0[:maxInd]),fractionP[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']  
 



if S==2:
    fractionP2=dataFile['fractionP2'][fieldGoodArray][sort]
    fractionM2=dataFile['fractionM2'][fieldGoodArray][sort]
    fractions=np.append(np.append(np.append(np.append(fractionM2[:maxInd],fractionM[:maxInd]),fraction0[:maxInd]),fractionP[:maxInd]),fractionP2[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']   
 
fractions=fractions.reshape(2*S+1,fraction0[:maxInd].size).transpose()   

figure=plt.figure()
pan=figure.add_subplot(1,1,1)
#pan.plot(tList,fractionP2+1,'co')
#pan.plot(tList,fractionP+0.5,'bo')
#pan.plot(tList,fractionM-0.5,'ro')
#pan.plot(tList,fractionM2-1,'mo')
pan.plot(tList,fraction0,'go')
pan.set_xlabel('Hold time [s]')
pan.set_ylabel('Fractions')



cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

kgrid=np.fft.fftfreq(tList.size,d=(tList[1]-tList[0]))
fft0=np.fft.fft(fraction0)
fftP=np.fft.fft(fractionP)
fftM=np.fft.fft(fractionM)
fftP2=np.fft.fft(fractionP2)
fftM2=np.fft.fft(fractionM2)

figure=plt.figure()
pan=figure.add_subplot(1,1,1)
pan.plot(kgrid,fft0*np.conjugate(fft0),cDict[0])
pan.plot(kgrid,fftP*np.conjugate(fftP),cDict[1])
pan.plot(kgrid,fftM*np.conjugate(fftM),cDict[-1])
pan.plot(kgrid,fftP2*np.conjugate(fftP2),cDict[2])
pan.plot(kgrid,fftM2*np.conjugate(fftM2),cDict[-2])

figure=plt.figure()
pan=figure.add_subplot(1,1,1)
pan.set_title(filename)
pan.plot(kgrid,fft0*np.conjugate(fft0)+fftP*np.conjugate(fftP)+fftM*np.conjugate(fftM)+fftP2*np.conjugate(fftP2)+fftM2*np.conjugate(fftM2),'k-')
pan.set_xlabel('Frequency [Hz]')
pan.set_ylabel('Sum of square moduli of FFt for all spin states')