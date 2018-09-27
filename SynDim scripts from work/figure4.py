# -*- coding: utf-8 -*-
"""
Created on Mon May 15 17:20:20 2017

@author: dng5
"""

import matplotlib.pyplot as plt
#import numpy as np



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np


#rcParams['axes.labelsize'] = 10
#rcParams['xtick.labelsize'] = 8
#rcParams['ytick.labelsize'] = 8
#rcParams['legend.fontsize'] = 8
#
#rcParams['pdf.fonttype'] = 42 # True type fonts
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Computer Modern Roman']
#rcParams['text.usetex'] = True
#rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
#
#rcParams['axes.linewidth'] = 0.75
#rcParams['lines.linewidth'] = 0.75
#
#rcParams['xtick.major.size'] = 3      # major tick size in points
#rcParams['xtick.minor.size'] = 2      # minor tick size in points
#rcParams['xtick.major.width'] = 0.75       # major tick width in points
#rcParams['xtick.minor.width'] = 0.75      # minor tick width in points
#
#rcParams['ytick.major.size'] = 3      # major tick size in points
#rcParams['ytick.minor.size'] = 2      # minor tick size in points
#rcParams['ytick.major.width'] = 0.75       # major tick width in points
#rcParams['ytick.minor.width'] = 0.75      # minor tick width in points

ColorMap = "afmhot"

datafileF1p=np.load('08Mar2017_F1_chern_1.npz')
datafileF1m=np.load('07Mar2017_F1_chern_-1.npz')
datafileF10=np.load('09Mar2017_Rf_Corrected.npz')
datafileF2p=np.load('28Feb2017_F2_chern_1.npz')
datafileF2m=np.load('27Feb2017_F2_chern_-1.npz')
datafileF20=np.load('22Mar2017_Rf_Corrected.npz')
width=4
lw=15
figure = plt.figure()
figure.clear()
figure.set_size_inches(6.5,3.5)

gs = gridspec.GridSpec(2,3)

gs.update(left=0.1, right=0.95, top=0.95, bottom=0.15)

qlistFit=datafileF1p['qlistFit']
fractionP=datafileF1p['fractionP']
fraction0=datafileF1p['fraction0']
fractionM=datafileF1p['fractionM']
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
pan2=figure.add_subplot(gs[0])
pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionP,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='o',linewidths=lw)
pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='o',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionM,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='o',linewidths=lw)
pan2.set_ylim([-3,3])
pan2.set_ylabel(r'Lattice site along $\hat{y}$')
pan2.xaxis.set_ticks([])

qlist=datafileF10['qlist']
fractionP=datafileF10['fractionP']
fraction0=datafileF10['fraction0']
fractionM=datafileF10['fractionM']

pan2=figure.add_subplot(gs[1])
pan2.scatter(qlist,[1 for j in range(qlist.size)],c=fractionP,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[0 for j in range(qlist.size)],c=fraction0,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[-1 for j in range(qlist.size)],c=fractionM,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.set_ylim([-3,3])
pan2.xaxis.set_ticks([])
pan2.yaxis.set_ticks([])

qlistFit=datafileF1m['qlistFit']
fractionP=datafileF1m['fractionP']
fraction0=datafileF1m['fraction0']
fractionM=datafileF1m['fractionM']
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
pan2=figure.add_subplot(gs[2])
pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionP,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionM,s=width,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=lw)
pan2.set_ylim([-3,3])
pan2.xaxis.set_ticks([])
pan2.yaxis.set_ticks([])

qlistFit=datafileF2p['qlistFit']
fractionP2=datafileF2p['fractionP2']
fractionP=datafileF2p['fractionP']
fraction0=datafileF2p['fraction0']
fractionM=datafileF2p['fractionM']
fractionM2=datafileF2p['fractionM2']
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
pan2=figure.add_subplot(gs[3])
pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionP,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionM,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[2 for j in range(qlistFit[inRange].size)],c=fractionP2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-2 for j in range(qlistFit[inRange].size)],c=fractionM2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
pan2.set_ylabel(r'Lattice site along $\hat{y}$')

qlist=datafileF20['qlist']
fractionP2=datafileF20['fractionP2']
fractionP=datafileF20['fractionP']
fraction0=datafileF20['fraction0']
fractionM=datafileF20['fractionM']
fractionM2=datafileF20['fractionM2']
pan2=figure.add_subplot(gs[4])
pan2.scatter(qlist,[1 for j in range(qlist.size)],c=fractionP,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[0 for j in range(qlist.size)],c=fraction0,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[-1 for j in range(qlist.size)],c=fractionM,s=wiFalsedth,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[2 for j in range(qlist.size)],c=fractionP2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlist,[-2 for j in range(qlist.size)],c=fractionM2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
pan2.yaxis.set_ticks([])

qlistFit=datafileF2m['qlistFit']
fractionP2=datafileF2m['fractionP2']
fractionP=datafileF2m['fractionP']
fraction0=datafileF2m['fraction0']
fractionM=datafileF2m['fractionM']
fractionM2=datafileF2m['fractionM2']
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
pan2=figure.add_subplot(gs[5])
pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionP,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionM,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[2 for j in range(qlistFit[inRange].size)],c=fractionP2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.scatter(qlistFit[inRange],[-2 for j in range(qlistFit[inRange].size)],c=fractionM2,s=width,vmin=0.0,vmax=0.5,cmap='Blues', marker='_',linewidths=lw)
pan2.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
pan2.yaxis.set_ticks([])




    
figure.show()
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure4.pdf', transparent=True)