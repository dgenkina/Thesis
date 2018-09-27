# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 13:22:15 2018

@author: dng5
"""
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize


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

cdictAll = {'red':  ((0.0, 1.0, 1.0),
                   (0.25, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.75, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),

         'blue':  ((0.0, 1.0, 1.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0, 1.0, 1.0))
        }
allColors = mpl.colors.LinearSegmentedColormap('allColors', cdictAll)
plt.register_cmap(cmap=allColors)

S=2
s=2*S+1
theoryFile=np.load('SynDimBandStructure_F2_n7_Chern1.npz')
kList = theoryFile['kList']
E = theoryFile['E']
m = theoryFile['m']
pops=theoryFile['pops']#[:,0,:]

modeTheory=np.zeros((kList.size,s))
sigmaModeTheory=np.zeros((kList.size,s))
fitGoodTheoryList=np.ones((kList.size,s),dtype=bool)
sigmaArrayTheory=np.ones(pops.shape)*0.0001
for bandInd in np.arange(s):
    for Kind in range(kList.size):
        try:
            modeTheory[Kind,bandInd],sigmaModeTheory[Kind,bandInd]=getMode(S,pops[Kind,bandInd,:],sigmaArrayTheory[Kind,bandInd,:])
        except RuntimeError:
            fitGoodTheoryList[Kind,bandInd]=False

figure=plt.figure(figsize=(3.73,2.4))
panel=figure.add_subplot(1,1,1)
plt.subplots_adjust(top=0.8,bottom=0.15,left=0.15)


for i in range(s):   
    d=panel.scatter(kList/2.0,E[:,i]/0.0775 + 6.0,c=modeTheory[:,i],cmap='allColors',vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-0.5,0.5)
 #   panel.tick_params(labelsize=20)
#panel.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
panel.set_ylabel(r'Energy [$t_x$]')
panel.set_ylim(-4.0,0.0)
panel.set_yticks([0.0,-1.0,-2.0,-3.0,-4.0])
panel.set_xticks([-0.5,-0.25,0,0.25,0.5])
panel.set_xticklabels([])
cbaxes = figure.add_axes([0.15, 0.82, 0.75, 0.05]) 
cbar = plt.colorbar(d,orientation='horizontal',ticks=np.arange(-S,S+1),cax = cbaxes)  

cbar.ax.tick_params(top=True,bottom=False,labelbottom=False, labeltop=True) 
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r'Modal position $\bar{m}$')#,labelpad=-28)
  
#plt.tight_layout()
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure2v7a.jpg',dpi=400)
np.savez('BandStructureFigureData',kList=kList,E=E,modeTheory=modeTheory)