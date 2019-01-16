# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 15:33:04 2017

@author: dng5
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams

rcParams['axes.labelsize'] = 10
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8

rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

rcParams['axes.linewidth'] = 0.75
rcParams['lines.linewidth'] = 0.75
rcParams['lines.markersize'] = 2

rcParams['xtick.major.size'] = 3      # major tick size in points
rcParams['xtick.minor.size'] = 2      # minor tick size in points
rcParams['xtick.major.width'] = 0.75       # major tick width in points
rcParams['xtick.minor.width'] = 0.75      # minor tick width in points

rcParams['ytick.major.size'] = 3      # major tick size in points
rcParams['ytick.minor.size'] = 2      # minor tick size in points
rcParams['ytick.major.width'] = 0.75       # major tick width in points
rcParams['ytick.minor.width'] = 0.75      # minor tick width in points


#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,7.0))
gs = gridspec.GridSpec(2,3)
gs.update(left=0.15, right=0.95, top=0.9, bottom = 0.1,wspace = 0.1,hspace=0.3)

#Read in positive kick, negative kick, no kick data
S=2
date='22Mar2017'
posFile=np.load('C:/Users/swooty/Documents/Thesis Data/Rf correction data/'+date+'_posKick.npz')
posMag=posFile['mag']
posNum=posFile['numPoints']
posSigmaMag=posFile['sigmaMag']
posFractionP=posFile['fractionP']
posSigmaP=posFile['sigmaP']/np.sqrt(posNum)
posFraction0=posFile['fraction0']
posSigma0=posFile['sigma0']/np.sqrt(posNum)
posFractionM=posFile['fractionM']
posSigmaM=posFile['sigmaM']/np.sqrt(posNum)

if S==2:
    posFractionP2=posFile['fractionP2']
    posSigmaP2=posFile['sigmaP2']/np.sqrt(posNum)
    posFractionM2=posFile['fractionM2']
    posSigmaM2=posFile['sigmaM2']/np.sqrt(posNum)
posTime=posFile['tlist']
posQ=posFile['qlistFit']
negFile=np.load('C:/Users/swooty/Documents/Thesis Data/Rf correction data/'+date+'_negKick.npz')
negNum=negFile['numPoints']
negMag=negFile['mag']
negSigmaMag=negFile['sigmaMag']

negFractionP=negFile['fractionP']
negSigmaP=negFile['sigmaP']/np.sqrt(negNum)
negFraction0=negFile['fraction0']
negSigma0=negFile['sigma0']/np.sqrt(negNum)
negFractionM=negFile['fractionM']
negSigmaM=negFile['sigmaM']/np.sqrt(negNum)

if S==2:
    negFractionP2=negFile['fractionP2']
    negSigmaP2=negFile['sigmaP2']/np.sqrt(negNum)
    negFractionM2=negFile['fractionM2']
    negSigmaM2=negFile['sigmaM2']/np.sqrt(negNum)
negTime=negFile['tlist']
negQ=negFile['qlistFit']
noFile=np.load('C:/Users/swooty/Documents/Thesis Data/Rf correction data/'+date+'_noKick.npz')
noMag=noFile['mag']
noSigmaMag=noFile['sigmaMag']
noFractionP=noFile['fractionP']
noFraction0=noFile['fraction0']
noFractionM=noFile['fractionM']
if S==2:
    noFractionP2=noFile['fractionP2']
    noFractionM2=noFile['fractionM2']
noTime=noFile['tlist']
noQ=noFile['qlist']

#Plot uncorrected data for all three

pan1=fig.add_subplot(gs[0,0])
pan1.plot(negTime,negFractionP2,'co',label='mF=+2')
pan1.plot(negTime,negFractionP,'bo',label='mF=+1')
pan1.plot(negTime,negFraction0,'go',label='mF=0')
pan1.plot(negTime,negFractionM,'ro',label='mF=-1')
pan1.plot(negTime,negFractionM2,'mo',label='mF=-2')
pan1.set_xlabel('Time  [s]')
pan1.set_ylabel('Fractional populations')
pan1.set_title('(a) negative force')
plt.legend(loc=2)

pan2=fig.add_subplot(gs[0,1])
pan2.plot(posTime,posFractionP,'bo',label='mF=+1')
pan2.plot(posTime,posFraction0,'go',label='mF=0')
pan2.plot(posTime,posFractionM,'ro',label='mF=-1')
if S==2:
    pan2.plot(posTime,posFractionP2,'co',label='mF=+2')
    pan2.plot(posTime,posFractionM2,'mo',label='mF=-2')
pan2.set_xlabel('Time  [s]')
pan2.set_yticklabels([])
pan2.set_title('(b) positive force')

#qlistFit=datafileF1p['qlistFit']
#
pan3=fig.add_subplot(gs[0,2])
#fractionP=datafileF1p['fractionP']
#fraction0=datafileF1p['fraction0']
#fractionM=datafileF1p['fractionM']
#inRange=((qlistFit<1.0)&(qlistFit>-1.0))
#pan2=figure.add_subplot(gs[0])
#pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionP,s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0,s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionM,s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan3.set_yticklabels([])
pan3.plot(noTime,noFractionP,'bo',label='mF=+1')
pan3.plot(noTime,noFraction0,'go',label='mF=0')
pan3.plot(noTime,noFractionM,'ro',label='mF=-1')
if S==2:
    pan3.plot(noTime,noFractionP2,'co',label='mF=+2')
    pan3.plot(noTime,noFractionM2,'mo',label='mF=-2')
pan3.set_xlabel('Time  [s]')
pan3.set_title('(c) no force')



# CALCULATE CORRECTION FACTORS BASED ON SHEER FRACTIONAL DIFFERENCE
noSort=np.argsort(noTime)
noDeltP=noFractionP-noFractionP[noSort][0]
noDelt0=noFraction0-noFraction0[noSort][0]
noDeltM=noFractionM-noFractionM[noSort][0]

if S==2:
    noDeltP2=noFractionP2-noFractionP2[noSort][0]
    noDeltM2=noFractionM2-noFractionM2[noSort][0]


##LOAD FITTED FOR DETUNING FRACTIONS FOR NO KICK DATA     
#noFileFit=np.load('C:/Users/swooty/Documents/Thesis Data/Rf correction data/'+date+'_noKick_detuningFit.npz')
#undetunedFracs=noFileFit['undetunedFracs']
#fittedFractions=noFileFit['fittedFractions']
#noTime2=noFileFit['tList']
#
##CALCULATE CORRECTION FACTORS BASED ON DETUNING
#noDeltP=fittedFractions[:,1+S]-undetunedFracs[1+S]
#noDelt0=fittedFractions[:,S]-undetunedFracs[S]
#noDeltM=fittedFractions[:,S-1]-undetunedFracs[S-1]
#
#if S==2:
#    noDeltP2=fittedFractions[:,2+S]-undetunedFracs[2+S]
#    noDeltM2=fittedFractions[:,S-2]-undetunedFracs[S-2]

#GET CORRECTED VERSIONS OF POSITIVE AND NEGATIVE KICK DATA, BASED ON TIME 
posFracPCorr=np.zeros(posMag.size)
posFrac0Corr=np.zeros(posMag.size)
posFracMCorr=np.zeros(posMag.size)
if S==2:
    posFracP2Corr=np.zeros(posMag.size)
    posFracM2Corr=np.zeros(posMag.size)
    
for i,t in enumerate(posTime):
    j=np.where(noTime==t)[0][0]
    posFracPCorr[i]=posFractionP[i]-noDeltP[j]
    posFrac0Corr[i]=posFraction0[i]-noDelt0[j]
    posFracMCorr[i]=posFractionM[i]-noDeltM[j]
    if S==2:
        posFracP2Corr[i]=posFractionP2[i]-noDeltP2[j]
        posFracM2Corr[i]=posFractionM2[i]-noDeltM2[j]
    
negFracPCorr=np.zeros(negMag.size)
negFrac0Corr=np.zeros(negMag.size)
negFracMCorr=np.zeros(negMag.size)
if S==2:
    negFracP2Corr=np.zeros(negMag.size)
    negFracM2Corr=np.zeros(negMag.size)
for i,t in enumerate(negTime):
    j=np.where(noTime==t)[0][0]
    negFracPCorr[i]=negFractionP[i]-noDeltP[j]
    negFrac0Corr[i]=negFraction0[i]-noDelt0[j]
    negFracMCorr[i]=negFractionM[i]-noDeltM[j]
    if S==2:
        negFracP2Corr[i]=negFractionP2[i]-noDeltP2[j]
        negFracM2Corr[i]=negFractionM2[i]-noDeltM2[j]

#PUT TOGETHER CORRECTED POSITIVE AND NEGATIVE FIT DATA 
qlist=np.append(posQ,negQ)
inRange=((qlist<1.0)&(qlist>-1.0))
if S==2:
    fractionP2=np.append(posFracP2Corr,negFracP2Corr)
    sigmaP2=np.append(posSigmaP2,negSigmaP2)
    fractionM2=np.append(posFracM2Corr,negFracM2Corr)
    sigmaM2=np.append(posSigmaM2,negSigmaM2)
fractionP=np.append(posFracPCorr,negFracPCorr)
sigmaP=np.append(posSigmaP,negSigmaP)
fraction0=np.append(posFrac0Corr,negFrac0Corr)
sigma0=np.append(posSigma0,negSigma0)
fractionM=np.append(posFracMCorr,negFracMCorr)
sigmaM=np.append(posSigmaM,negSigmaM)

#PLOT CORRECTED DATA AS A FUNCTION OF CRYSTAL MOMENTUM

pan=fig.add_subplot(gs[1,2])
pan.errorbar(qlist[inRange], fractionP[inRange],yerr=sigmaP[inRange],fmt='bo',label=r'$m_F=+1$')
pan.errorbar(qlist[inRange], fraction0[inRange],yerr=sigma0[inRange],fmt='go',label=r'$m_F=0$')
pan.errorbar(qlist[inRange], fractionM[inRange],yerr=sigmaM[inRange],fmt='ro',label=r'$m_F=-1$')
if S==2:
    pan.errorbar(qlist[inRange], fractionP2[inRange],yerr=sigmaP2[inRange],fmt='co',label=r'$m_F=+2$')
    pan.errorbar(qlist[inRange], fractionM2[inRange],yerr=sigmaM2[inRange],fmt='mo',label=r'$m_F=-2$')
pan.set_yticklabels([])
pan.set_xlabel(r'crystal momentum [$k_L$]')
pan.set_title('(f) combined')
#legend()

#fig2=plt.figure()
#pan2=fig2.add_subplot(1,1,1)
#pan2.scatter(qlist[inRange],[1 for j in range(qlist[inRange].size)],c=fractionP[inRange],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#pan2.scatter(qlist[inRange],[0 for j in range(qlist[inRange].size)],c=fraction0[inRange],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#pan2.scatter(qlist[inRange],[-1 for j in range(qlist[inRange].size)],c=fractionM[inRange],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#if S==2:
#    pan2.scatter(qlist[inRange],[2 for j in range(qlist[inRange].size)],c=fractionP2[inRange],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#    pan2.scatter(qlist[inRange],[-2 for j in range(qlist[inRange].size)],c=fractionM2[inRange],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
#pan2.set_ylabel('Synthetic lattice site')
#pan2.set_xlabel(r'Crystal momentum [$k_L$]')

#PLOT CORRECTED DATA AS A FUNCTION OF TIME, INDIVIDUALLY FOR DIFFERENT KICKS

pan1=fig.add_subplot(gs[1,0])
pan1.plot(negTime,negFracPCorr,'bo',label='mF=+1')
pan1.plot(negTime,negFrac0Corr,'go',label='mF=0')
pan1.plot(negTime,negFracMCorr,'ro',label='mF=-1')
if S==2:
    pan1.plot(negTime,negFracP2Corr,'co',label='mF=+2')
    pan1.plot(negTime,negFracM2Corr,'mo',label='mF=-2')
pan1.set_xlabel('Time  [s]')
pan1.set_ylabel('Corrected fractional populations')
pan1.set_title('(d) negative force')

pan2=fig.add_subplot(gs[1,1])
pan2.plot(posTime,posFracPCorr,'bo',label='mF=+1')
pan2.plot(posTime,posFrac0Corr,'go',label='mF=0')
pan2.plot(posTime,posFracMCorr,'ro',label='mF=-1')
if S==2:
    pan2.plot(posTime,posFracP2Corr,'co',label='mF=+2')
    pan2.plot(posTime,posFracM2Corr,'mo',label='mF=-2')
pan2.set_xlabel('Time  [s]')
pan2.set_title('(e) positive force')
pan2.set_yticklabels([])


plt.savefig('RfCorrection.pdf',transparent=True)

##PLOT UNCORRECTED MAGNETIZATIONS AS A FUNCTION OF TIME FOR ALL THREE KICKS
#figure=plt.figure()
#pan=figure.add_subplot(1,1,1)
#pan.errorbar(posTime,posMag,yerr=posSigmaMag,fmt='bo',label='positive kick')
#pan.errorbar(negTime,negMag,yerr=negSigmaMag,fmt='go',label='negative kick')
#pan.errorbar(noTime,noMag,yerr=noSigmaMag,fmt='ro',label='no kick')
#pan.set_xlabel('Time [s]')
#pan.set_ylabel('Magnetization')
#plt.legend()
#
##CORRECT MAGNETIZATION BY SUBTRACTING NO KICK MAGNETIZATION AND PLOT
#posMagCorr=np.zeros(posMag.size)
#posSigmaCorr=np.zeros(posSigmaMag.size)
#for i,t in enumerate(posTime):
#    j=np.where(noTime==t)[0][0]
#    posMagCorr[i]=posMag[i]-noMag[j]
#    posSigmaCorr=np.sqrt(posSigmaMag[i]**2.0+noSigmaMag[j]**2.0)
#    
#figure=plt.figure()
#pan=figure.add_subplot(1,1,1)
#pan.errorbar(posTime,posMagCorr,yerr=posSigmaCorr,fmt='bo')
#pan.set_xlabel('Time [s]')
#pan.set_ylabel('Corrected magnetization, positive kick')
#
#
#negMagCorr=np.zeros(negMag.size)
#negSigmaCorr=np.zeros(negSigmaMag.size)
#for i,t in enumerate(negTime):
#    j=np.where(noTime==t)[0][0]
#    negMagCorr[i]=negMag[i]-noMag[j]
#    negSigmaCorr=np.sqrt(negSigmaMag[i]**2.0+noSigmaMag[j]**2.0)
#    
#figure=plt.figure()
#pan=figure.add_subplot(1,1,1)
#pan.errorbar(negTime,negMagCorr,yerr=negSigmaCorr,fmt='bo')
#pan.set_xlabel('Time [s]')
#pan.set_ylabel('Corrected magnetization, negative kick')
#
#Filename=date+'_Rf_Corrected'
#np.savez(Filename,qlist=qlist[inRange],fractionM=fractionM[inRange],fraction0=fraction0[inRange],fractionP=fractionP[inRange],
#         sigmaM=sigmaM[inRange],sigma0=sigma0[inRange],sigmaP=sigmaP[inRange])
#         
#if S==2:
#    np.savez(Filename,qlist=qlist[inRange],fractionM2=fractionM2[inRange],fractionM=fractionM[inRange],fraction0=fraction0[inRange],
#             fractionP=fractionP[inRange],fractionP2=fractionP2[inRange],sigmaM2=sigmaM2[inRange],
#             sigmaM=sigmaM[inRange],sigma0=sigma0[inRange],sigmaP=sigmaP[inRange],sigmaP2=sigmaP2[inRange])