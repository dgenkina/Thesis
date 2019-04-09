# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 17:29:40 2019

@author: swooty
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import matplotlib as mpl
import lineFit2
import figureChenOfOmegaSize

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

rcParams['xtick.major.size'] = 3      # major tick size in points
rcParams['xtick.minor.size'] = 2      # minor tick size in points
rcParams['xtick.major.width'] = 0.75       # major tick width in points
rcParams['xtick.minor.width'] = 0.75      # minor tick width in points

rcParams['ytick.major.size'] = 3      # major tick size in points
rcParams['ytick.minor.size'] = 2      # minor tick size in points
rcParams['ytick.major.width'] = 0.75       # major tick width in points
rcParams['ytick.minor.width'] = 0.75      # minor tick width in points

cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
k_L = 2.0*np.pi/lambdaL

plotBool = False

filestart = 'TDSEkickF2experparams'
numlist = np.arange(6.0,16.0)
forcelist = (lambdaL/Erecoil)*hbar*2.0*k_L/(numlist*1e-3)
filelist = [filestart+str(np.round(numlist[ind],1))+'.npz' for ind in range(numlist.size)]
filelist = np.array(filelist)

fig = plt.figure()
fig.clear()
fig.set_size_inches(6.5,3.5)
gs = gridspec.GridSpec(1,2)
gs.update(left=0.15, right=0.95, top=0.85, bottom = 0.2, wspace=0.25, hspace=0.25)
pan = fig.add_subplot(gs[0,0])

chernLine=np.zeros(filelist.size)
delChernLine=np.zeros(filelist.size)
chernDio=np.zeros(filelist.size)
delChernDio=np.zeros(filelist.size)
for ind,filename in enumerate(filelist):
    TDSEfile = np.load(filename)
    qlistTDSE = (TDSEfile['tgrid']*hbar/Erecoil-TDSEfile['latticeRampOnt']-TDSEfile['ramanRampOnt'])/(TDSEfile['tgrid'][-1]*hbar/Erecoil-TDSEfile['latticeRampOnt']-TDSEfile['ramanRampOnt'])
    modeTDSE = TDSEfile['modeList'][qlistTDSE>0]
    fracs = TDSEfile['fracs'][qlistTDSE>0]
    qlistTDSE = qlistTDSE[qlistTDSE>0]
    qlistTDSE = np.append(-qlistTDSE[::-1], qlistTDSE)
    modeTDSE = np.append(-modeTDSE[::-1], modeTDSE)
    fracs = np.append(fracs[::-1,::-1],fracs, axis=0)
    
    pan.plot(qlistTDSE/2.0,modeTDSE,label=np.round(forcelist[ind],3))
    
    A,B,dA,dB = lineFit2.lineFit(qlistTDSE,modeTDSE,'q theory','m theory',plot=plotBool, bounds=(np.array([-np.inf,-0.001]),np.array([np.inf,0.001])))
    chernLine[ind] = A*2.0/3.0
    delChernLine[ind] = dA*2.0/3.0
    
    indM=qlistTDSE.size/6
    ind0=qlistTDSE.size/2
    indP=5*qlistTDSE.size/6
    
    xPth,dxPth=figureChenOfOmegaSize.getMaxFitTheory(qlistTDSE,fracs[:,1],np.ones(fracs[:,1].size)*0.001,indC=indP,plot=plotBool)
    x0th,dx0th=figureChenOfOmegaSize.getMaxFitTheory(qlistTDSE,fracs[:,2],np.ones(fracs[:,2].size)*0.001,indC=ind0,plot=plotBool)
    xMth,dxMth=figureChenOfOmegaSize.getMaxFitTheory(qlistTDSE,fracs[:,3],np.ones(fracs[:,3].size)*0.001,indC=indM,plot=plotBool)
    Ath,Bth,dAth,dBth=lineFit2.lineFit(np.array([xMth,x0th,xPth]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxMth,dx0th,dxPth]),plot=plotBool)
    chernDio[ind] = Ath*2.0/3.0
    delChernDio[ind] = dAth*2.0/3.0
pan.set_xlabel(r'Crystal momentum $q_x$ [2$k_L$]')
pan.set_ylabel(r'Modal position $\bar{m}$')
#plt.legend()

pan=fig.add_subplot(gs[0,1])
pan.errorbar(forcelist,chernLine,yerr=delChernLine,fmt='b-')
pan.errorbar(forcelist,chernDio,yerr=delChernDio,fmt='g-')
pan.set_xlabel(r'Force $F_x$ [$E_L/\lambda_L$]')
pan.set_ylabel(r'Inferred Chern number')
plt.savefig('C:\Users\swooty\Documents\Work\Bloch Osc Paper\Figures\ChernForForceRangeFig.pdf', transparent=True)

