# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 11:58:01 2019

@author: swooty
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
import SynDimBandStructureGeneral as sdbs

Flat = True
omega=0.5
delta=0.0
epsilon = 0.0
U = 4.4
n = 7
S = 2
m0 = 0
c=1064.0/790.0

fileroot = 'C:/Users/swooty/Documents/Thesis Data/Final Bloch Osc data/'
filename10 = 'SynDimBandStructure_F7_n13_Chern1Flat.npz'
filename2 = 'SynDimBandStructure_F2_n7_Chern1Flat.npz'
file10=np.load(fileroot+filename10)
file2=np.load(fileroot+filename2)
#kList2,E2,pops2, m2 = sdbs.plotSynDimBandStructGen(omega, delta, epsilon, U, 7,2, m0,c=c, Flat=Flat)
#
#kList20,E20,pops20, m20 = sdbs.plotSynDimBandStructGen(omega, delta, epsilon, U, 15,20, m0,c=c, Flat=Flat)

plt.close(1)
fig = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,9)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.2,hspace=0.2)

S = 2
s = 2*S+1
panel=fig.add_subplot(gs[0,0:4])
for i in range(s):   
    d=panel.scatter(file2['kList'],file2['E'][:,i],c=file2['m'][:,i],vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-1.0,1.0)

S = 7
s = 2*S+1
panel=fig.add_subplot(gs[0,4:8])
for i in range(s):   
    d=panel.scatter(file10['kList'],file10['E'][:,i],c=file10['m'][:,i],vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-1.0,1.0)

#make colorbar for F=2
cbar_axes = fig.add_subplot(gs[0,8])
cbar = fig.colorbar(d,cax=cbar_axes)
cbar_axes.yaxis.set_ticks_position('left')