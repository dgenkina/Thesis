# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 17:14:19 2017

@author: dng5
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import lineFit2
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8


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
    
filenamelist=['07Mar2017_F1_chern_-1.npz','27Feb2017_F2_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz','22Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz','28Feb2017_F2_chern_1.npz']
filenamelist2=['22Mar2017_Rf_Corrected.npz',
               '28Feb2017_F2_chern_1.npz',
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
clist2=[0,1,-1]
titleList=['(a)','(d)','(b)','(e)','(c)','(f)']
titleList2=['(a)','(b)','(c)']
fig = plt.figure()
fig.clear()
fig.set_size_inches(3.5,3.5)
gs = gridspec.GridSpec(3,2)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15, hspace=0.5, wspace=0.5)
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
    
    
       
    pan=fig.add_subplot(gs[i*2])
    pan.errorbar(qlist[fitGoodList]/2.0, mode[fitGoodList],yerr=sigmaMode[fitGoodList],fmt='o',color='Grey',markersize=2)
    pan.plot(kList[fitGoodTheoryList]/2.0,modeTheory[fitGoodTheoryList],'r-',linewidth=1)
    
    
    qrangeList = np.linspace(0.2,1.0,20)
    Alist = np.zeros(qrangeList.size)
    dAlist = np.zeros(qrangeList.size)
    numPointsList = np.zeros(qrangeList.size)

    for ind,qrange in enumerate(qrangeList):
        qrangeBool=((qlist[fitGoodList]<qrange) & (qlist[fitGoodList]>-qrange))
        Alist[ind],B,dAlist[ind],dB = lineFit2.lineFit(qlist[fitGoodList][qrangeBool], mode[fitGoodList][qrangeBool],'q data','m data',yerr=sigmaMode[fitGoodList][qrangeBool],plot=False,errorbars=True, bounds=(np.array([-np.inf,-0.001]),np.array([np.inf,0.001])))
        numPointsList[ind] = qlist[fitGoodList][qrangeBool].size
    #pan.set_xlabel(r'$q_x$ [$k_L$]',size=14)
    #pan.set_ylabel(r'$\langle m \rangle$',size=14)
    qrangeListTh = np.linspace(0.01,1.0,25)
    AlistTh = np.zeros(qrangeListTh.size)
    dAlistTh = np.zeros(qrangeListTh.size)
    numPointsListTh = np.zeros(qrangeListTh.size)
    for ind,qrange in enumerate(qrangeListTh):
        krangeBool=((kList[fitGoodTheoryList]<qrange) & (kList[fitGoodTheoryList]>-qrange))
        AlistTh[ind],B,dAlistTh[ind],dB = lineFit2.lineFit(kList[fitGoodTheoryList][krangeBool],modeTheory[fitGoodTheoryList][krangeBool],'q theory','m theory',plot=False, bounds=(np.array([-np.inf,-0.001]),np.array([np.inf,0.001])))
        numPointsListTh[ind] = kList[fitGoodTheoryList][krangeBool].size
    #if clist2[i]!=0:
#    axins = inset_axes(pan,
#                       width="20%", # width = 30% of parent_bbox
#                       height="30%", # height : 1 inch
#                       bbox_to_anchor=bbox[i])
    pan2=fig.add_subplot(gs[i*2+1])
    pan2.errorbar(qrangeList/2.0,Alist*2.0/3.0,yerr=dAlist*2.0/3.0,fmt='o',color='Grey',markersize=2)
    pan2.errorbar(qrangeListTh/2.0,AlistTh*2.0/3.0,yerr=dAlistTh*2.0/3.0,fmt='r-')
   # axins.set_ylabel('Slope')
   # axins.set_yticks([0.5*clist2[i],1.0*clist2[i]])
   # axins.set_xlabel(r'Max. $|q|$ in fit range' )
    pan2.set_xticks([0.0, 0.25, 0.5])
    pan2.set_xticklabels([])
    pan2.set_xlim(0.0,0.5)
 #   pan2.yaxis.tick_right()
    pan2.set_yticks([-0.5+0.5*clist2[i],0.0+0.5*clist2[i],0.5+0.5*clist2[i]])
    pan2.set_ylim(-0.5+0.5*clist2[i],0.5+0.5*clist2[i])
    if i==1:
        #pan2.yaxis.set_label_position("right")
        pan2.set_ylabel(r'Conductivity ${\tilde\sigma}_{\rm{H}}$')
    if i==2:

        pan2.set_xticklabels(['0', '0.25', '0.5'])
        pan2.set_xlabel(r'Max $|q_x|$ in fit')
#        
    extraFig=plt.figure()
    extraPan=extraFig.add_subplot(111)
    extraPan.errorbar(qrangeList,Alist*2.0/3.0,yerr=dAlist/np.sqrt(numPointsList),fmt='k-',label='experiment')
    extraPan.errorbar(qrangeListTh,AlistTh*2.0/3.0,yerr=dAlistTh/np.sqrt(numPointsListTh),fmt='r-',label='theory')
    extraPan.set_ylabel('Slope')
   # extraPan.set_yticks([0.2*clist2[i],0.6*clist2[i],1.0*clist2[i]])
    extraPan.set_xlabel(r'Max. $|q|$ in fit range' )
    extraPan.set_xticks([0.2,0.6,1.0])
    extraPan.set_xlim(0.2,1.0)
    plt.legend()
    plt.savefig('C:/Users/swooty/Documents/Thesis/SynDim scripts from work/slopeOfRange%i.png' %(i),dpi=400)
#    
    pan.plot(kList/2.0,lineFit2.line(kList,3.0*clist2[i]/2.0,0),'k--',linewidth=2)
    
    dataMode[str(clist2[i])]=mode[fitGoodList]
    dataModeError[str(clist2[i])]=sigmaMode[fitGoodList]
    theoryMode[str(clist2[i])]=modeTheory[fitGoodTheoryList]
    qlistData[str(clist2[i])]=qlist[fitGoodList]/2.0
    qlistTheory[str(clist2[i])]=kList[fitGoodTheoryList]/2.0
    
    pan.set_ylim(-S,S)
    pan.set_xlim(-0.5,0.5)
    pan.set_xticks([])
    pan.set_yticks(np.arange(-S,S+1))
    
    pan.set_xticks([-0.5,-0.25,0,0.25,0.5])
    pan.set_xticklabels([])
    #    pan.set_title(titleList2[i],loc='left')
    #    pan.set_ylabel(r'$\bar{m}$')
    if i==2:
        pan.set_xticklabels([-0.5,-0.25,0,0.25,0.5])
        pan.set_xlabel(r'Crystal momentum $q_x$ [$2k_L$]')
    if i==1:
        pan.set_ylabel(r'Modal position $\bar{m}$')
    
   
fig.show()
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/magnetizationAll3.pdf',dpi=400)
#plt.savefig('Z:/My Documents/papers etc/talks and posters/magnetizationAll3.jpg',dpi=400)
#totNum=dataFile['numM2'][fieldGoodArray]+dataFile['numM'][fieldGoodArray]+dataFile['num0'][fieldGoodArray]+dataFile['numP'][fieldGoodArray]+dataFile['numP2'][fieldGoodArray]
#time=dataFile['tlist'][fieldGoodArray]
qlist=dataFile['qlist'][fieldGoodArray]
#bgndAvg=dataFile['bgndAvg'][fieldGoodArray]

numPcorr=fractionP#/((1.5e-5)*yCenters**2.0-0.00538*yCenters+1.329)
num0corr=fraction0#/((1.6e-5)*yCenters**2.0-0.00583*yCenters+1.296)
numMcorr=fractionM#/((1.4e-5)*yCenters**2.0-0.00473*yCenters+1.266)

#np.savez('magnetizationPlotData',dataMode0=dataMode['0'],dataModeP=dataMode['1'],dataModeM=dataMode['-1'],
#         dataModeError0=dataModeError['0'],dataModeErrorP=dataModeError['1'],dataModeErrorM=dataModeError['-1'],
#         theoryMode0=theoryMode['0'],theoryModeP=theoryMode['1'],theoryModeM=theoryMode['-1'],
#         qlistData0=qlistData['0'],qlistDataP=qlistData['1'],qlistDataM=qlistData['-1'],
#         qlistTheory0=qlistTheory['0'],qlistTheoryP=qlistTheory['1'],qlistTheoryM=qlistTheory['-1'])
