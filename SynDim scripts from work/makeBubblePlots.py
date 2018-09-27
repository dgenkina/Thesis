# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:16:20 2017

@author: dng5
"""
import matplotlib.pyplot as plt
#import numpy as np



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np

datafileF1p=np.load('08Mar2017_F1_chern_1.npz')
datafileF1m=np.load('07Mar2017_F1_chern_-1.npz')
datafileF10=np.load('09Mar2017_Rf_Corrected.npz')
datafileF2p=np.load('28Feb2017_F2_chern_1.npz')
datafileF2m=np.load('27Feb2017_F2_chern_-1.npz')
datafileF20=np.load('22Mar2017_Rf_Corrected.npz')

qlistFit=datafileF2m['qlist']
fractionP=datafileF2m['fractionP']
fraction0=datafileF2m['fraction0']
fractionM=datafileF2m['fractionM']
S=2

if S==2:
    fractionP2=datafileF2m['fractionP2']
    fractionM2=datafileF2m['fractionM2']
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
qlist=qlistFit[inRange]
Half=np.remainder(np.arange(qlist.size),2)==0
figure=plt.figure()
figure.set_size_inches(6.5,3.5)
pan2=figure.add_subplot(111)
pan2.scatter(qlist[Half],[1 for j in range(qlist[Half].size)],c='b',edgecolors='b', marker='o',linewidths=np.sqrt(fractionP[Half])*10.0)
pan2.scatter(qlist[Half],[0 for j in range(qlist[Half].size)],c='g',edgecolors='g',marker='o',linewidths=np.sqrt(fraction0[Half])*10.0)
pan2.scatter(qlist[Half],[-1 for j in range(qlist[Half].size)],c='r',edgecolors='r',marker='o',linewidths=np.sqrt(fractionM[Half])*10.0)
if S==2:
    pan2.scatter(qlist[Half],[2 for j in range(qlist[Half].size)],c='c',edgecolors='c', marker='o',linewidths=np.sqrt(fractionP2[Half])*10.0)
    pan2.scatter(qlist[Half],[-2 for j in range(qlist[Half].size)],c='m',edgecolors='m',marker='o',linewidths=np.sqrt(fractionM2[Half])*10.0)
    
pan2.set_xlabel(r'Crystal momentum [$k_L$]',size=24)
pan2.set_ylabel(r'Lattice site m', size=24)
pan2.tick_params(labelsize=20)

  
plt.tight_layout()