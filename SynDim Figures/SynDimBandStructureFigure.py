# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 10:24:50 2018

@author: swooty
"""

import numpy as np
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from matplotlib import cm
import SynDimBandStructureGeneral as sd

from matplotlib import rc, rcParams


# Use me to produce pdf files
# from matplotlib.backends.backend_pgf import FigureCanvasPgf
# matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

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

#make figure
plt.close(1)
fig2 = plt.figure(1,figsize=(6.2,5.0))
gs = gridspec.GridSpec(2,10)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,hspace = 0.3,wspace=4.5)

#Get F=1 Omega=0 band structure
c=1064.0/790.0
S=1
kListF1O0,EF1O0,popsF1O0, mF1O0 = sd.plotSynDimBandStructGen(0.0, 0.0, 
                                                             0.02, 5.0, 7,S,0,c=c,
                           kList=np.linspace(-2.0,2.0,600),plot=False)

#make F=1 Omega=0 panel
panel=fig2.add_subplot(gs[0,0:5])
for i in range(2*S+1):   
    d=panel.scatter(kListF1O0,EF1O0[:,i],c=mF1O0[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-2.0,2.0)
panel.set_ylabel(r'Energy [$E_L$]')
panel.set_title('(a)')
panel.axvline(x=1.0,c='k')
panel.axvline(x=-1.0,c='k')

#Get F=1 Omega=0.5 band structure
c=1064.0/790.0
S=1
kListF1O05,EF1O05,popsF1O05, mF1O05 = sd.plotSynDimBandStructGen(0.5, 0.0, 
                                                             0.02, 5.0, 7,S,0,c=c,
                           kList=np.linspace(-1.0,1.0,600),plot=False)

#make F=1 Omega=0.5 panel
panel=fig2.add_subplot(gs[0,5:9])
for i in range(2*S+1):   
    d=panel.scatter(kListF1O05,EF1O05[:,i],c=mF1O05[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-1.0,1.0)

panel.set_title('(b)')

#make colorbar for F=1
cbar_axes = fig2.add_subplot(gs[0,9])
cbar = fig2.colorbar(d,cax=cbar_axes,ticks=np.array([-1,0,1]))
cbar_axes.yaxis.set_ticks_position('left')

#Get F=2 Omega=0 band structure
c=1064.0/790.0
S=2
kListF2O0,EF2O0,popsF2O0, mF2O0 = sd.plotSynDimBandStructGen(0.0, 0.0, 
                                                             0.02, 5.0, 7,S,0,c=c,
                           kList=np.linspace(-2.0,2.0,600),plot=False)

#make F=2 Omega=0 panel
panel=fig2.add_subplot(gs[1,0:5])
for i in range(2*S+1):   
    d=panel.scatter(kListF2O0,EF2O0[:,i],c=mF2O0[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)
panel.axvline(x=1.0,c='k')
panel.axvline(x=-1.0,c='k')
panel.set_xlim(-2,2)
panel.set_ylabel(r'Energy [$E_L$]')
panel.set_xlabel(r'$q$ [$k_L$]')
panel.set_title('(c)')

#Get F=2 Omega=0.5 band structure
c=1064.0/790.0
S=2
kListF2O05,EF2O05,popsF2O05, mF2O05 = sd.plotSynDimBandStructGen(0.5, 0.0, 
                                                             0.02, 5.0, 7,S,0,c=c,
                           kList=np.linspace(-1.0,1.0,600),plot=False)

#make F=2 Omega=0 panel
panel=fig2.add_subplot(gs[1,5:9])
for i in range(2*S+1):   
    d=panel.scatter(kListF2O05,EF2O05[:,i],c=mF2O05[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-1.0,1.0)
panel.set_xlabel(r'$q$ [$k_L$]')
panel.set_title('(d)')

#make colorbar for F=2
cbar_axes = fig2.add_subplot(gs[1,9])
cbar = fig2.colorbar(d,cax=cbar_axes,ticks=np.array([-2,-1,0,1,2]))
cbar_axes.yaxis.set_ticks_position('left')

plt.savefig('SynDimBandStructure.pdf', transparent=True)
