# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:32:12 2017

@author: dng5
"""
import numpy as np
import matplotlib.pyplot as plt


dataFile=np.load('23Jun2017_files_94-273_sw2_0.01.npz')

S=2

checkCoh=False
imbal=dataFile['imbalArray']
signalGood=dataFile['signalGood']
cutoff=0.2
fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
fileGood=(fieldGoodArray)
if checkCoh:
    coherent=dataFile['coherentList']
    fileGood=(fieldGoodArray & coherent)
times=dataFile['tlist'][fileGood]
sort=np.argsort(times)

qlist=dataFile['qlist'][fileGood]

figure=plt.figure()
pan=figure.add_subplot(1,1,1)
pan.plot(times,qlist,'bo')
pan.set_xlabel('pulse time [s]')
pan.set_ylabel('Quasimomentum')

fracP2=dataFile['fractionP2'][fileGood]
fracP=dataFile['fractionP'][fileGood]
frac0=dataFile['fraction0'][fileGood]
fracM=dataFile['fractionM'][fileGood]
fracM2=dataFile['fractionM2'][fileGood]

numBins=50
bins=np.linspace(np.min(times),np.max(times),numBins+1)

fractionP=np.zeros(numBins)
fraction0=np.zeros(numBins)
fractionM=np.zeros(numBins)

sigmaP=np.zeros(numBins)
sigma0=np.zeros(numBins)
sigmaM=np.zeros(numBins)

time=np.zeros(numBins)
q=np.zeros(numBins)
sigmaQ=np.zeros(numBins)
if S==2:
    fractionP2=np.zeros(numBins)
    fractionM2=np.zeros(numBins)
    sigmaP2=np.zeros(numBins)
    sigmaM2=np.zeros(numBins)
for ind in range(numBins):
    binInds=np.where((times>=bins[ind]) & (times<bins[ind+1]))[0]
    
    fractionP[ind]=np.average(fracP[binInds])
    sigmaP[ind]=np.std(fracP[binInds])/np.sqrt(binInds.size)
    fraction0[ind]=np.average(frac0[binInds])
    sigma0[ind]=np.std(frac0[binInds])/np.sqrt(binInds.size)
    fractionM[ind]=np.average(fracM[binInds])
    sigmaM[ind]=np.std(fracM[binInds])/np.sqrt(binInds.size)
    
    time[ind]=np.average(times[binInds])
    q[ind]=np.average(qlist[binInds])
    sigmaQ[ind]=np.std(qlist[binInds])/np.sqrt(binInds.size)
    if S==2:
        fractionP2[ind]=np.average(fracP2[binInds])
        sigmaP2[ind]=np.std(fracP2[binInds])/np.sqrt(binInds.size)
        fractionM2[ind]=np.average(fracM2[binInds])
        sigmaM2[ind]=np.std(fracM2[binInds])/np.sqrt(binInds.size)
    
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(time,fractionP,yerr=sigmaP, fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(time,fraction0,yerr=sigma0, fmt='go', label=r'$m_F$=0')
pan1.errorbar(time,fractionM,yerr=sigmaM, fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan1.errorbar(time,fractionP2,yerr=sigmaP2, fmt='co', label=r'$m_F$=+2')
    pan1.errorbar(time,fractionM2,yerr=sigmaM2, fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan1.set_xlabel('Time [s]')
pan1.set_ylabel('Fractional populations')
#legend()
    
    
    
    