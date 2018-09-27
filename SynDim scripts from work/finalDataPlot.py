# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 17:19:45 2017

@author: dng5
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import matplotlib as mpl
import lineFit2
import figureChenOfOmegaSize
import scipy.ndimage as nd


rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['axes.titlesize'] = 8


rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['text.latex.preamble'] = [r'\usepackage{cmbright}',r'\usepackage{amsmath}']


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

cdictC = {'red':   [(0.0,  1.0, 1.0), 
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.5, 0.5)],

         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  0.5, 0.5)]}
                   
cyans = mpl.colors.LinearSegmentedColormap('Cyans', cdictC)
plt.register_cmap(cmap=cyans)        

cdictM = {'red':   [(0.0,  1.0, 1.0), 
                   (1.0,  0.5, 0.5)],

         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,   0.5,  0.5)]}
                   
magentas = mpl.colors.LinearSegmentedColormap('Magentas', cdictM)
plt.register_cmap(cmap=magentas)

cmapDict={}
cmapDict[-2]='Magentas'
cmapDict[-1]='Reds'
cmapDict[0]='Greens'
cmapDict[1]='Blues'
cmapDict[2]='Cyans'



def parabola(xlist,A,B,x0):
    return A*(xlist-x0)**2.0+B
    
def gaussian(xlist,sigma,A,x0,B):
    return A*np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))/(sigma*np.sqrt(2.0*np.pi))+B
    
def fitParabola(xlist,ylist,yerr,plot=False):
    (A,B,x0), cov=optimize.curve_fit(parabola,xlist,ylist,p0=(-4.0,np.max(ylist),xlist[np.where(ylist==np.max(ylist))[0][0]]),sigma=yerr,absolute_sigma=True)
    (dA,dB,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        xforplot=np.linspace(np.min(xlist),np.max(xlist),num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.plot(xforplot,parabola(xforplot,A,B,x0))
        pan.set_title('A=%.2f +/-%.3f , B=%.2f +/-%.3f , x0=%.2f +/-%.3f ' %(A,dA,B,dB,x0,dx0))
        
    return x0,dx0
def fitGaussian(xlist,ylist,yerr,plot=False):
    (sigma,A,x0,B), cov=optimize.curve_fit(gaussian,xlist,ylist,p0=(3.0,np.max(ylist),xlist[np.where(ylist==np.max(ylist))[0][0]],np.min(ylist)))#,sigma=yerr,absolute_sigma=False)
    (dSigma,dA,dx0,dB)=np.sqrt(np.diag(cov))
    
    if plot:
        xforplot=np.linspace(np.min(xlist),np.max(xlist),num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.plot(xforplot,gaussian(xforplot,sigma,A,x0,B))
        pan.set_title(r'$\sigma$=%.2f+/-%.3f, A=%.2f+/-%.3f, x0=%.2f+/-%.3f, B=%.2f+/-%.3f' %(sigma,dSigma,A,dA,x0,dx0,B,dB))
        
    return x0,dx0
    
def getMaxFit(qlist,fracList,sigma,plot=False):
    sort=np.argsort(qlist)
    indP=np.where(fracList[sort]==np.max(fracList[sort]))[0][0]
    indPmin=indP-20
    indPmax=indP+20
    if indP < 20:
        indPmin=0
    if indP > fractionP.size-20:
        indPmax = fractionP.size
        
    
    xP,dxP=fitParabolags = gridspec.GridSpec(3,1, height_ratios=[15, 5,3]) 
    return xP,dxP
    
def getMaxFitTheory(qlist,fracList,sigma,plot=False):
    w=15
    sort=np.argsort(qlist)
    indP=np.where(fracList[sort]==np.max(fracList[sort]))[0][0]
    indCent=int(qlist.size/2)
    shift=indCent-indP
    qlistShift=np.roll(qlist[sort],shift)
    if shift>0:
        qlistShift[0:shift]=qlistShift[0:shift]-2.0
    elif shift<0:
        qlistShift[shift:qlist.size]=qlistShift[shift:qlist.size]+2.0
    qlistShift=qlistShift[indCent-w:indCent+w+1]
    fracShift=np.roll(fracList[sort],shift)[indCent-w:indCent+w+1]
    sigmaShift=np.roll(sigma[sort],shift)[indCent-w:indCent+w+1]
    
    xP,dxP=fitParabola(qlistShift,fracShift,sigmaShift,plot=plot)
    if xP>1.0:
        xP=xP-2.0
    if xP<-1.0:
        xP=xP+2.0
        
    return xP,dxP

filenamelist=['08Mar2017_F1_chern_1.npz']    
filenamelist=['07Mar2017_F1_chern_-1.npz','27Feb2017_F2_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz','22Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz','28Feb2017_F2_chern_1.npz']
plotBool=True
width=15
lw=10
Slist=[1]
Slist=[1,2,1,2,1,2]
clist=[1]
clist=[-1,-1,0,0,1,1]
titleList=['(a)','(d)','(b)','(e)','(c)','(f)']
mainFig = plt.figure()
mainFig.clear()
mainFig.set_size_inches(3.5,4.0)
gs = gridspec.GridSpec(3,2)
gs.update(left=0.15, right=0.95, top=0.85, bottom = 0.1, wspace=0.15, hspace=0.25)

filename='27Feb2017_F2_chern_-1.npz'
i=1

for i,filename in enumerate(filenamelist):
    print filename
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
    
    sort=np.argsort(qlist)
    
    if clist[i] != 0:
        if clist[i]>0:
            indM=kList.size/6
            ind0=kList.size/2
            indP=5*kList.size/6
            indMd=qlist.size/6
            ind0d=qlist.size/2
            indPd=5*qlist.size/6
        elif clist[i]<0:
            indM=5*kList.size/6
            ind0=kList.size/2
            indP=kList.size/6
            indMd=5*qlist.size/6
            ind0d=qlist.size/2
            indPd=qlist.size/6
        else:
            indM=kList.size/2
            ind0=kList.size/2
            indP=kList.size/2
            indMd=qlist.size/2
            ind0d=qlist.size/2
            indPd=qlist.size/2
        if S==1:
            xMth,dxMth=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,0],np.ones(pops[:,0].size)*0.001,indC=indM,plot=plotBool)
            x0th,dx0th=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,1],np.ones(pops[:,1].size)*0.001,indC=ind0,plot=plotBool)
            xPth,dxPth=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,2],np.ones(pops[:,2].size)*0.001,indC=indP,plot=plotBool)
        elif S==2:
            xMth,dxMth=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,1],np.ones(pops[:,1].size)*0.001,indC=indM,plot=plotBool)
            x0th,dx0th=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,2],np.ones(pops[:,2].size)*0.001,indC=ind0,plot=plotBool)
            xPth,dxPth=figureChenOfOmegaSize.getMaxFitTheory(kList,pops[:,3],np.ones(pops[:,3].size)*0.001,indC=indP,plot=plotBool)
        xP,dxP=figureChenOfOmegaSize.getMaxFitTheory(qlist,fractionP,sigmaP,indC=indPd,plot=plotBool)
        x0,dx0=figureChenOfOmegaSize.getMaxFitTheory(qlist,fraction0,sigma0,indC=ind0d,plot=plotBool)
        xM,dxM=figureChenOfOmegaSize.getMaxFitTheory(qlist,fractionM,sigmaM,indC=indMd,plot=plotBool)
    
        print 'S = ' + str(S)
        A,B,dA,dB=lineFit2.lineFit(np.array([xM,x0,xP]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxM,dx0,dxP]),plot=plotBool)
        print 'chern = '+ str(A*2/3) + '+/-'+str(dA*2/3)
           
        Ath,Bth,dAth,dBth=lineFit2.lineFit(np.array([xMth,x0th,xPth]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxMth,dx0th,dxPth]),plot=plotBool)
        print 'chern theory = '+ str(Ath*2/3) + '+/-'+str(dAth*2/3)
    
    fracs=np.array([fractionM,fraction0,fractionP])
    if S==2:
        fracs=np.array([fractionM2,fractionM,fraction0,fractionP,fractionP2])
    pan=mainFig.add_subplot(gs[i])
    for ind in range(2*S+1):
        pan.scatter(qlist/2,[ind-S for j in range(qlist.size)],c=fracs[ind],s=width,vmin=0.0,vmax=1.0,cmap=cmapDict[ind-S], marker='_',linewidths=lw)
    if clist[i] != 0:
        pan.errorbar(np.array([xM,x0,xP])/2,np.array([-1,0,1]),yerr=np.array([dxM,dx0,dxP]),fmt='ko',markersize=4)
        pan.plot(qlist/2,A*qlist+B,'k-')
        
    pan.xaxis.set_ticks([])
    pan.yaxis.set_ticks([])

    pan.set_ylim(-3,3)
    pan.set_xlim(-0.5,0.5)
    pan.set_yticks(np.arange(-2,3))
    pan.set_yticklabels([])
    pan.set_xticks([-0.5,-0.25,0,0.25,0.5])
    pan.set_xticklabels([])

    if ((i==4) or (i==5)):
        pan.set_xlabel(r'Crystal momentum $q_x$ [$2k_L$]')
        pan.set_xticklabels([-0.5,-0.25,0,0.25,0.5])

#        
    if (np.remainder(i,2)==0):
        pan.set_ylabel(r'site $m$')
        pan.set_yticklabels([-2,-1,0,1,2])


mainFig.show()
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/barcodeWithSlope.pdf', dpi=500)

