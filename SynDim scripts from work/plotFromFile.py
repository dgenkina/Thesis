# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:31:29 2016

@author: dng5
"""
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize

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

rcParams['xtick.major.size'] = 3      # major tick size in points
rcParams['xtick.minor.size'] = 2      # minor tick size in points
rcParams['xtick.major.width'] = 0.75       # major tick width in points
rcParams['xtick.minor.width'] = 0.75      # minor tick width in points

rcParams['ytick.major.size'] = 3      # major tick size in points
rcParams['ytick.minor.size'] = 2      # minor tick size in points
rcParams['ytick.major.width'] = 0.75       # major tick width in points
rcParams['ytick.minor.width'] = 0.75      # minor tick width in points

def gaussian(xlist,sigma,x0):
    return np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))/(sigma*np.sqrt(2.0*np.pi))

def getMode(S,fracList,sigmaList,plot=False):
    mFlist=np.arange(2*S+1)-S
    mF=np.lib.pad(mFlist, (1,1), 'reflect', reflect_type='odd')
    frac=np.lib.pad(fracList,(1,1),'constant',constant_values=(0,0))
    yerr=np.lib.pad(sigmaList,(1,1),'constant',constant_values=(0.0001,0.0001))
    
    (sigma,x0), cov=optimize.curve_fit(gaussian,mF,frac,p0=(S/1.0,float(mF[np.where(np.max(frac)==frac)[0][0]])),sigma=yerr,absolute_sigma=False)
    (dSigma,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        mFforplot=np.linspace(mF[0],mF[-1],num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(mF,frac,yerr=yerr,fmt='bo')
        pan.plot(mFforplot,gaussian(mFforplot,sigma,x0))
    return x0, dx0
    
filenamelist=['07Mar2017_F1_chern_-1.npz','27Feb2017_F2_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz','22Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz','28Feb2017_F2_chern_1.npz']
#filename='07Mar2017_F1_chern_-1.npz'
kick='pos'
filelist=np.arange(354,384)
xlabel = 'oscillation time [s]'
#S=1
Slist=[1,2,1,2,1,2]
clist=[-1,-1,0,0,1,1]
fig = plt.figure()
fig.clear()
fig.set_size_inches(3.0,5.0)
gs = gridspec.GridSpec(3,2)
gs.update(left=0.18, right=0.95, top=0.95, bottom = 0.15)

for i,filename in enumerate(filenamelist):
    dataFile=np.load(filename)
    S=Slist[i]
    #imbal=dataFile['imbalArray']
    #signalGood=dataFile['signalGood']
    cutoff=0.2
    fieldGoodArray=np.ones(dataFile['qlist'].size,dtype=bool)#fieldGoodArray=((np.abs(imbal)<cutoff) & signalGood)#
    
    qlist=dataFile['qlist']
    fractionP=dataFile['fractionP'][fieldGoodArray]
    fraction0=dataFile['fraction0'][fieldGoodArray]
    fractionM=dataFile['fractionM'][fieldGoodArray]
    sigmaP=dataFile['sigmaP'][fieldGoodArray]
    sigma0=dataFile['sigma0'][fieldGoodArray]
    sigmaM=dataFile['sigmaM'][fieldGoodArray]
    mag=fractionP-fractionM
    sigmaMag=np.sqrt(sigmaP**2.0+sigmaM**2.0)
    fracArray=np.array((fractionM,fraction0,fractionP))
    sigmaArray=np.array((sigmaM,sigma0,sigmaP))
    
    if S==2:
        fractionP2=dataFile['fractionP2'][fieldGoodArray]
        fractionM2=dataFile['fractionM2'][fieldGoodArray]
        sigmaP2=dataFile['sigmaP2'][fieldGoodArray]
        sigmaM2=dataFile['sigmaM2'][fieldGoodArray]
        fracArray=np.array((fractionM2,fractionM,fraction0,fractionP,fractionP2))
        sigmaArray=np.array((sigmaM2,sigmaM,sigma0,sigmaP,sigmaP2))
        mag=2.0*fractionP2+fractionP-fractionM-2.0*fractionM2
        sigmaMag=np.sqrt(4.0*sigmaP2**2.0+sigmaP**2.0+sigmaM**2.0+4.0*sigmaM2**2.0)
    theoryFile=np.load('SynDimBandStructure_F'+str(S)+'_n7_Chern'+str(clist[i])+'.npz')
    kList=theoryFile['kList']
    pops=theoryFile['pops'][:,0,:]
        
    mode=np.zeros(fraction0.size)
    sigmaMode=np.zeros(fraction0.size)
    fitGoodList=np.ones(fraction0.size,dtype=bool)
    for ind in range(fraction0.size):
        try:
            mode[ind],sigmaMode[ind]=getMode(S,fracArray[:,ind],sigmaArray[:,ind])
        except RuntimeError:
            fitGoodList[ind]=False
        if sigmaMode[ind]>0.5:
            fitGoodList[ind]=False
    
    modeTheory=np.zeros(kList.size)
    sigmaModeTheory=np.zeros(kList.size)
    fitGoodTheoryList=np.ones(kList.size,dtype=bool)
    sigmaArrayTheory=np.ones(pops.shape)*0.0001
    for ind in range(kList.size):
        try:
            modeTheory[ind],sigmaModeTheory[ind]=getMode(S,pops[ind,:],sigmaArrayTheory[ind,:])
        except RuntimeError:
            fitGoodTheoryList[ind]=False
       
    pan=fig.add_subplot(gs[i])
    pan.errorbar(qlist[fitGoodList], mode[fitGoodList],yerr=sigmaMode[fitGoodList],fmt='bo')
    #pan.set_xlabel(r'$q_x$ [$k_L$]',size=14)
    #pan.set_ylabel(r'$\langle m \rangle$',size=14)
    pan.plot(kList[fitGoodTheoryList],modeTheory[fitGoodTheoryList],'b-')
    pan.tick_params(labelsize=12)
    pan.set_ylim(-S,S)
    pan.set_xlim(-1.0,1.0)
    pan.set_xticks([-1.0,0,1.0])
    pan.set_yticks([-S,0,S])
    if np.remainder(i,2)==0:
        pan.set_ylabel(r'Expected site along $\hat{s}$')
    if ((i==4) or (i==5)):
        pan.set_xlabel(r'$q_x$ [$k_L$]')
  #  plt.tight_layout()

fig.show()
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/magnetizationAll.pdf', transparent=True)

#totNum=dataFile['numM2'][fieldGoodArray]+dataFile['numM'][fieldGoodArray]+dataFile['num0'][fieldGoodArray]+dataFile['numP'][fieldGoodArray]+dataFile['numP2'][fieldGoodArray]
#time=dataFile['tlist'][fieldGoodArray]
qlist=dataFile['qlist'][fieldGoodArray]
#bgndAvg=dataFile['bgndAvg'][fieldGoodArray]

numPcorr=fractionP#/((1.5e-5)*yCenters**2.0-0.00538*yCenters+1.329)
num0corr=fraction0#/((1.6e-5)*yCenters**2.0-0.00583*yCenters+1.296)
numMcorr=fractionM#/((1.4e-5)*yCenters**2.0-0.00473*yCenters+1.266)
#
#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(time,totNum,'co')
#pan1.set_xlabel(xlabel)
#pan1.set_ylabel('total counted od')
#
#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(time,dataFile['num0'][fieldGoodArray],'co')
#pan1.set_xlabel(xlabel)
#pan1.set_ylabel('total counted od in mF=0')


#sort=np.argsort(time)
#qlist=qlist[sort]*2.0
#filelist=filelist[fieldGoodArray][sort]
#if kick=='neg':
#    for i in range(qlist.size-1):
#        if qlist[i+1]<qlist[i]-1.0:
#            qlist[i+1]=qlist[i+1]+2.0
#            
#if kick=='pos':
#    for i in range(qlist.size-1):
#        if qlist[i+1]>qlist[i]+0.8:
#            qlist[i+1]=qlist[i+1]-2.0
#numP2=dataFile['numP2'][fieldGoodArray]
#numP=dataFile['numP'][fieldGoodArray]
#num0=dataFile['num0'][fieldGoodArray]
#numM=dataFile['numM'][fieldGoodArray]
#numM2=dataFile['numM2'][fieldGoodArray]

#A,B,dA,dB=lineFit(time[sort],qlist,'pulse time',r'quasimomentum [$k_L$]')
#qlist=qlist-B


#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(time[sort],qlist,'co')
#pan1.set_xlabel(xlabel)
#pan1.set_ylabel('quasimomentum [k_L]')
#
#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(bgndAvg,fractionP2,'co')
#pan1.plot(bgndAvg,fractionP,'bo')
#pan1.plot(bgndAvg,fraction0,'go')
#pan1.plot(bgndAvg,fractionM,'ro')
#pan1.plot(bgndAvg,fractionM2,'mo')
#pan1.set_xlabel('average background counts')
#pan1.set_ylabel('Fractional pops')
#pan1.set_title(filename)
#legend()

#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(time,fractionP2,'co')
#pan1.plot(time,fractionP,'bo')
#pan1.plot(time,fraction0,'go')
#pan1.plot(time,fractionM,'ro')
#pan1.plot(time,fractionM2,'mo')
#pan1.set_xlabel(r'Time [s]')
#pan1.set_ylabel('Fractional population')
#pan1.set_title(filename)
#legend()
#
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.plot(qlist,fractionP2,'co',label=r'$m_F=2')
pan1.plot(qlist,fractionP,'bo',label=r'$m_F=1')
pan1.plot(qlist,fraction0,'go',label=r'$m_F=0')
pan1.plot(qlist,fractionM,'ro',label=r'$m_F=-1')
pan1.plot(qlist,fractionM2,'mo',label=r'$m_F=-2')
pan1.set_xlabel(r'Quasimomentum [$k_L$]')
pan1.set_ylabel('Fractional population')
pan1.set_title(filename)
pan1.set_ylim(0,1.2)
pan1.set_xlim(-1,1)
legend()