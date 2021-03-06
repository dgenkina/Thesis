# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 19:21:12 2019

@author: swooty
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import lineFit2
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rcParams['axes.labelsize'] = 10
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10


rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage[cm]{sfmath}']

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
rcParams['xtick.direction'] = 'in' #ticks inside graph
rcParams['ytick.direction'] = 'in' #ticks inside graph

def gaussian(xlist,sigma,x0):
    return np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))/(sigma*np.sqrt(2.0*np.pi))

def getMode(S,fracList,sigmaList,centGuess='nan',plot=False):
    pad=1
    mFlist=np.arange(2*S+1)-S
    mF=np.lib.pad(mFlist, (pad,pad), 'reflect', reflect_type='odd')
    frac=np.lib.pad(fracList, (pad,pad),'constant',constant_values=(0,0))
    yerr=np.lib.pad(sigmaList, (pad,pad),'constant',constant_values=(0.02,0.02))
    
    if centGuess=='nan':
        centGuess=float(mF[np.where(np.max(frac)==frac)[0][0]])
    
    (sigma,x0), cov=optimize.curve_fit(gaussian,mF,frac,p0=(S/(S+1.0),centGuess),sigma=yerr,absolute_sigma=False)
    (dSigma,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        print centGuess
        mFforplot=np.linspace(mF[0],mF[-1],num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(mF,frac,yerr=yerr,fmt='bo')
        pan.plot(mFforplot,gaussian(mFforplot,sigma,x0))
    return x0, dx0


hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
    
filenamelist=['07Mar2017_F1_chern_-1.npz','27Feb2017_F2_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz','22Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz','28Feb2017_F2_chern_1.npz']
filenamelist2=['28Feb2017_F2_chern_1.npz',
               '22Mar2017_Rf_Corrected.npz',
              '27Feb2017_F2_chern_-1.npz']      
filenamelist1=['07Mar2017_F1_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz'] 
filename='27Feb2017_F2_chern_-1.npz'
i=2
kick='pos'
filelist=np.arange(354,384)
xlabel = 'oscillation time [s]'

Slist=[1,2,1,2,1,2]
clist=[-1,-1,0,0,1,1]
clist2=[1,0,-1]
titleList=['(a)','(d)','(b)','(e)','(c)','(f)']
titleList2=['(a)','(b)','(c)']
fig = plt.figure()
fig.clear()
fig.set_size_inches(3.5, 3.0)
gs = gridspec.GridSpec(3,1)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15, hspace=0.3, wspace=0.5)
dataMode={}
dataModeError={}
theoryMode={}
qlistData={}
qlistTheory={}

bbox={}
bbox[0] = [-180,215,400,100]
bbox[1] = [35,76,200,50]#[-180,215,400,100]
bbox[2] = [-92,15,200,50]#[-180,35,400,100]
fileroot = 'C:/Users/swooty/Documents/Thesis Data/Final Bloch Osc data/'
TDSEfile = np.load('TDSEkickF2experparams.npz')
qlistTDSE = (TDSEfile['tgrid']*hbar/Erecoil-TDSEfile['latticeRampOnt']-TDSEfile['ramanRampOnt'])/(TDSEfile['tgrid'][-1]*hbar/Erecoil-TDSEfile['latticeRampOnt']-TDSEfile['ramanRampOnt'])
modeTDSE = TDSEfile['modeList'][qlistTDSE>0]
qlistTDSE = qlistTDSE[qlistTDSE>0]
qlistTDSE = np.append(-qlistTDSE[::-1], qlistTDSE)
modeTDSE = np.append(-modeTDSE[::-1], modeTDSE)

for i,filename in enumerate(filenamelist2):
    dataFile=np.load(fileroot + filename)
    S=2#Slist[i]
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
    theoryFile=np.load(fileroot+'SynDimBandStructure_F'+str(S)+'_n7_Chern'+str(clist2[i])+'.npz')
    kList=theoryFile['kList']
    pops=theoryFile['pops'][:,0,:]
    
    modeTheory=np.zeros(kList.size)
    sigmaModeTheory=np.zeros(kList.size)
    fitGoodTheoryList=np.ones(kList.size,dtype=bool)
    sigmaArrayTheory=np.ones(pops.shape)*0.0001
    for ind in range(kList.size):
        try:
            modeTheory[ind],sigmaModeTheory[ind]=getMode(S,pops[ind,:],sigmaArrayTheory[ind,:])
        except RuntimeError:
            fitGoodTheoryList[ind]=False
    
        
    mode=np.zeros(fraction0.size)
    sigmaMode=np.zeros(fraction0.size)
    fitGoodList=np.ones(fraction0.size,dtype=bool)
    for ind in range(fraction0.size):
        centTheory = modeTheory[np.where(np.min(np.abs(kList-qlist[ind]))==np.abs(kList-qlist[ind]))[0][0]]
        centTheory = np.round(centTheory, decimals=3)
        try:
            mode[ind],sigmaMode[ind]=getMode(S,fracArray[:,ind],sigmaArray[:,ind], centGuess=centTheory)
        except RuntimeError:
            fitGoodList[ind]=False
        if sigmaMode[ind]>0.5:
            fitGoodList[ind]=False
            
    indP=np.where(np.min(np.abs(qlist-2/3))==np.abs(qlist-2/3))[0][0]
    ind0=np.where(np.min(np.abs(qlist))==np.abs(qlist))[0][0]
    indM=np.where(np.min(np.abs(qlist+2/3))==np.abs(qlist+2/3))[0][0]
    chern=((mode[indP]-mode[ind0])+(mode[ind0]-mode[indM]))/2.0
    sigmaChern=((sigmaMode[indP]-sigmaMode[ind0])+(sigmaMode[ind0]-sigmaMode[indM]))/2.0
    
    #print 'Chern measured = ' +str(chern)+'+/-'+str(sigmaChern)
    
    
       
    pan=fig.add_subplot(gs[i])
    #data points
    pan.errorbar(qlist[fitGoodList], mode[fitGoodList],yerr=sigmaMode[fitGoodList],fmt='o',color='Grey',markersize=4)
    #Lowest band theory lines
   # pan.plot(kList[fitGoodTheoryList],modeTheory[fitGoodTheoryList],'r-',linewidth=1)
    
    
    qrangeList = np.linspace(0.2,1.0,20)
    Alist = np.zeros(qrangeList.size)
    dAlist = np.zeros(qrangeList.size)
    numPointsList = np.zeros(qrangeList.size)

    for ind,qrange in enumerate(qrangeList):
        qrangeBool=((qlist[fitGoodList]<qrange) & (qlist[fitGoodList]>-qrange))
        Alist[ind],B,dAlist[ind],dB = lineFit2.lineFit(qlist[fitGoodList][qrangeBool], mode[fitGoodList][qrangeBool],'q data','m data',yerr=sigmaMode[fitGoodList][qrangeBool],plot=False,errorbars=True, bounds=(np.array([-np.inf,-0.001]),np.array([np.inf,0.001])))
        numPointsList[ind] = qlist[fitGoodList][qrangeBool].size

    qrangeListTh = np.linspace(0.01,1.0,25)
    AlistTh = np.zeros(qrangeListTh.size)
    dAlistTh = np.zeros(qrangeListTh.size)
    numPointsListTh = np.zeros(qrangeListTh.size)
    for ind,qrange in enumerate(qrangeListTh):
        krangeBool=((kList[fitGoodTheoryList]<qrange) & (kList[fitGoodTheoryList]>-qrange))
        AlistTh[ind],B,dAlistTh[ind],dB = lineFit2.lineFit(kList[fitGoodTheoryList][krangeBool],modeTheory[fitGoodTheoryList][krangeBool],'q theory','m theory',plot=False, bounds=(np.array([-np.inf,-0.001]),np.array([np.inf,0.001])))
        numPointsListTh[ind] = kList[fitGoodTheoryList][krangeBool].size

    #linear IQHE-like theory
    pan.plot(kList,lineFit2.line(kList,3.0*clist2[i]/2.0,0),'k--',linewidth=2)
    #full TDSE theory
   # pan.plot(qlistTDSE,modeTDSE*clist2[i],'b-',linewidth=2)
    
    dataMode[str(clist2[i])]=mode[fitGoodList]
    dataModeError[str(clist2[i])]=sigmaMode[fitGoodList]
    theoryMode[str(clist2[i])]=modeTheory[fitGoodTheoryList]
    qlistData[str(clist2[i])]=qlist[fitGoodList]/2.0
    qlistTheory[str(clist2[i])]=kList[fitGoodTheoryList]/2.0
    
    pan.set_ylim(-S,S)
    pan.set_xlim(-1.0,1.0)
    pan.set_xticks([])
    pan.set_yticks(np.arange(-S,S+1))
    
    pan.set_xticks([-1.0,-0.5,0,0.5,1.0])
    pan.set_xticklabels([])
    if i==2:
        pan.set_xticklabels([-1.0,-0.5,0,0.5,1.0])
        pan.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
    if i==1:
        pan.set_ylabel(r'Modal position $\bar{m}$')
    
   
fig.show()
plt.savefig('magnetizationAll1.png',dpi=400)

