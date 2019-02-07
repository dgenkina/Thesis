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
    
def gaussian(xlist,sigma,A,x0,B):
    return A*np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))+B

def fitGauss(xlist,ylist,yerr,plot=False):
    (sigma,A,x0,B), cov=optimize.curve_fit(gaussian,xlist,ylist,p0=(0.56,np.max(ylist),xlist[np.where(ylist==np.max(ylist))[0][0]],0.0),sigma=yerr,absolute_sigma=False)
    (dsigma,dA,dx0,dB)=np.sqrt(np.diag(cov))
    
    xforplot=np.linspace(np.min(xlist),np.max(xlist),num=100)
    yforplot = gaussian(xforplot,sigma,A,x0,B)
    
    if plot:
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.plot(xforplot,yforplot)
        pan.set_title('sigma=%.2f +/-%.3f , A=%.2f +/-%.3f , x0=%.2f +/-%.3f, B=%.2f +/-%.3f ' %(sigma,dsigma,A,dA,x0,dx0,B,dB))
       
    return x0,dx0,xforplot,yforplot
    
def getMaxFitTheory(qlist,fracList,sigma,indC='nan',plot=False):
    w=np.int(1.0*qlist.size/4.0)
    sort=np.argsort(qlist)
    if indC == 'nan':
        indC=np.where(fracList[sort]==np.max(fracList[sort]))[0][0]
    indCent=int(qlist.size/2)
    shift=indCent-indC
    qlistShift=np.roll(qlist[sort],shift)
    if shift>0:
        qlistShift[0:shift]=qlistShift[0:shift]-2.0
    elif shift<0:
        qlistShift[shift:qlist.size]=qlistShift[shift:qlist.size]+2.0
    qlistShift=qlistShift[indCent-w+1:indCent+w]

    fracShift=np.roll(fracList[sort],shift)[indCent-w+1:indCent+w]
    sigmaShift=np.roll(sigma[sort],shift)[indCent-w+1:indCent+w]
    
    xP,dxP,xforplot,yforplot=fitGauss(qlistShift,fracShift,sigmaShift,plot=plot)
 #   xP,dxP=getExpectationValue(qlistShift,fracShift,sigmaShift,plot=plot)
    if xP>1.0:
        xP=xP-2.0
    if xP<-1.0:
        xP=xP+2.0
        
    return xP,dxP,xforplot,yforplot
def line(x,A,B):
    return A*x+B
    
def lineFit(xvars,yvars,xlabel,ylabel,errorbars=False,yerr=0.0,plot=True,**kwargs):
    if errorbars:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,sigma=yerr,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        if plot:
            figure=plt.figure()
            pan=figure.add_subplot(1,1,1)
            pan.errorbar(xvars,yvars,yerr=yerr,fmt='bo')
            pan.plot(xrangefit,data_fitted,'b-')
            pan.set_title('Fit params in Ax+B, A='+str(np.round(A,3))+'+/-'+str(np.round(dA,4))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
            pan.set_xlabel(xlabel)
            pan.set_ylabel(ylabel)
    else:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        if plot:
            figure=plt.figure()
            pan=figure.add_subplot(1,1,1)
            pan.plot(xvars,yvars,'bo')
            pan.plot(xrangefit,data_fitted,'b-')
            pan.set_title('Fit params in Ax+B, A='+str(np.round(A,3))+'+/-'+str(np.round(dA,4))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
            pan.set_xlabel(xlabel)
            pan.set_ylabel(ylabel)
    return A,B,dA,dB,xrangefit,data_fitted

fileroot='C:/Users/swooty/Documents/Thesis Data/Final Bloch Osc data/'
filenamelist=[fileroot+'08Mar2017_F1_chern_1.npz']    
#filenamelist=[fileroot+'07Mar2017_F1_chern_-1.npz',fileroot+'27Feb2017_F2_chern_-1.npz',
#              fileroot+'09Mar2017_Rf_Corrected.npz',fileroot+'22Mar2017_Rf_Corrected.npz',
#              fileroot+'08Mar2017_F1_chern_1.npz',fileroot+'28Feb2017_F2_chern_1.npz']
plotBool=False
width=15
lw=10
Slist=[1]
#Slist=[1,2,1,2,1,2]
clist=[1]
#clist=[-1,-1,0,0,1,1]
filename=fileroot+'07Mar2017_F1_chern_-1.npz'#'08Mar2017_F1_chern_1.npz'
i=1

print filename
dataFile=np.load(filename)
S=1
c=-1
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
theoryFile=np.load(fileroot+'SynDimBandStructure_F'+str(S)+'_n7_Chern'+str(c)+'.npz')
kList=theoryFile['kList']
pops=theoryFile['pops'][:,0,:]

sort=np.argsort(qlist)

if c>0:
    indM=kList.size/6
    ind0=kList.size/2
    indP=5*kList.size/6
    indMd=qlist.size/6
    ind0d=qlist.size/2
    indPd=5*qlist.size/6

if c<0:
    indM=5*kList.size/6
    ind0=kList.size/2
    indP=kList.size/6
    indMd=5*qlist.size/6
    ind0d=qlist.size/2
    indPd=qlist.size/6

    
plt.close(1)
fig = plt.figure(1,figsize=(6.2,6.0))
gs = gridspec.GridSpec(2,2)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.2,hspace=0.2)
if S==1:
    xMth,dxMth,xfitMth,yfitMth=getMaxFitTheory(kList,pops[:,0],np.ones(pops[:,0].size)*0.001,indC=indM,plot=plotBool)
    x0th,dx0th,xfit0th,yfit0th=getMaxFitTheory(kList,pops[:,1],np.ones(pops[:,1].size)*0.001,indC=ind0,plot=plotBool)
    xPth,dxPth,xfitPth,yfitPth=getMaxFitTheory(kList,pops[:,2],np.ones(pops[:,2].size)*0.001,indC=indP,plot=plotBool)
elif S==2:
    xMth,dxMth,xfitMth,yfitMth=getMaxFitTheory(kList,pops[:,1],np.ones(pops[:,1].size)*0.001,indC=indM,plot=plotBool)
    x0th,dx0th,xfit0th,yfit0th=getMaxFitTheory(kList,pops[:,2],np.ones(pops[:,2].size)*0.001,indC=ind0,plot=plotBool)
    xPth,dxPth,xfitPth,yfitpth=getMaxFitTheory(kList,pops[:,3],np.ones(pops[:,3].size)*0.001,indC=indP,plot=plotBool)
xP,dxP,xfitP,yfitP=getMaxFitTheory(qlist,fractionP,sigmaP,indC=indPd,plot=plotBool)
x0,dx0,xfit0,yfit0=getMaxFitTheory(qlist,fraction0,sigma0,indC=ind0d,plot=plotBool)
xM,dxM,xfitM,yfitM=getMaxFitTheory(qlist,fractionM,sigmaM,indC=indMd,plot=plotBool)

pan = fig.add_subplot(gs[0,0])
pan.errorbar(qlist,fractionP,yerr=sigmaP,fmt='bo')
pan.plot(xfitP,yfitP,'b-')
pan.errorbar(qlist,fraction0,yerr=sigma0,fmt='go')
pan.plot(xfit0,yfit0,'g-')
pan.errorbar(qlist,fractionM,yerr=sigmaM,fmt='ro')
pan.plot(xfitM,yfitM,'r-')
pan.set_xlim(-1.25,1.25)
pan.set_ylabel('Fractional populations')
#pan.set_title('(a)')
pan.set_xticklabels([])

pan = fig.add_subplot(gs[0,1])
pan.plot(kList,pops[:,0],'ro')
pan.plot(xfitMth,yfitMth,'r-')
pan.plot(kList,pops[:,1],'go')
pan.plot(xfit0th,yfit0th,'g-')
pan.plot(kList,pops[:,2],'bo')
pan.plot(xfitPth,yfitPth,'b-')
pan.set_xlim(-1.25,1.25)
#pan.set_title('(b)')
pan.set_xticklabels([])

print 'S = ' + str(S)
A,B,dA,dB,xLine,yLine=lineFit(np.array([xM,x0,xP]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxM,dx0,dxP]),plot=plotBool)
print 'chern = '+ str(A*2/3) + '+/-'+str(dA*2/3)
   
Ath,Bth,dAth,dBth,xLineTh,yLineTh=lineFit(np.array([xMth,x0th,xPth]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxMth,dx0th,dxPth]),plot=plotBool)
print 'chern theory = '+ str(Ath*2/3) + '+/-'+str(dAth*2/3)
    

pan=fig.add_subplot(gs[1,0])
pan.errorbar(np.array([xM,x0,xP]),np.array([-1,0,1]),yerr=np.array([dxM,dx0,dxP]),fmt='ko')
pan.plot(xLine,yLine,'b-')
pan.set_xlim(-1.25,1.25)
#pan.set_title('Fit params in Ax+B, A='+str(np.round(A,3))+'+/-'+str(np.round(dA,4))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
pan.text(-0.5,-1.0,r'$C=$'+str(np.round(A*2/3,2)) + r'$\pm$'+str(np.round(dA*2/3,3)))
pan.set_yticks([-1,0,1])
pan.set_xlabel(r'crystal momentum [$k_L$]')
pan.set_ylabel(r'site $m$')
#pan.set_title('(c)')

pan=fig.add_subplot(gs[1,1])
pan.errorbar(np.array([xMth,x0th,xPth]),np.array([-1,0,1]),yerr=np.array([dxMth,dx0th,dxPth]),fmt='ko')
pan.plot(xLineTh,yLineTh,'b-')
pan.set_xlim(-1.25,1.25)
pan.set_yticks([-1,0,1])
pan.text(-0.5,-1.0,r'$C=$'+str(np.round(Ath*2/3,3)) + r'$\pm$'+str(np.round(dAth*2/3,4)))
#pan.set_title('Fit params in Ax+B, A='+str(np.round(Ath,3))+'+/-'+str(np.round(dAth,4))+', B='+str(np.round(Bth,3))+'+/-'+str(np.round(dBth,4)))
pan.set_xlabel(r'crystal momentum [$k_L$]')
#pan.set_title('(d)')

#plt.savefig('findingChern.pdf',transparent=True)