# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 10:24:50 2018

@author: swooty
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import SynDimBandStructureGeneral as sd

from matplotlib import rcParams


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
fig2 = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,10)
gs.update(left=0.13, right=0.95, top=0.9, bottom = 0.15,wspace=4.7)


#Get F=2 Omega=0 syn. dim. epsilon band structure
c=1064.0/790.0
S=2
kListF2bad,EF2bad,popsF2bad, mF2bad = sd.plotSynDimBandStructGen(0.0, 0.0, 
                                                             0.05, 4.4, 7,S,0,c=c,
                           kList=np.linspace(-1.0,1.0,600),plot=False)

#make F=2 Omega=0 syn. dim. epsilon panel
panel=fig2.add_subplot(gs[0,0:4])
for i in range(2*S+1):   
    d=panel.scatter(kListF2bad,EF2bad[:,i],c=mF2bad[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)

panel.set_xlim(-1,1)
panel.set_ylabel(r'Energy [$E_L$]')
panel.set_xlabel(r'$q$ [$k_L$]')
panel.set_title('(a)')

#Get F=2 Omega=0.5 bloch osc. epsilon band structure
c=1064.0/790.0
S=2
kListF2good,EF2good,popsF2good, mF2good = sd.plotSynDimBandStructGen(0.0, 0.0, 
                                                             0.02, 4.4, 7,S,0,c=c,
                           kList=np.linspace(-1.0,1.0,600),plot=False)

#make F=2 Omega=0 bloch osc. epsilon panel
panel=fig2.add_subplot(gs[0,4:8])
for i in range(2*S+1):   
    d=panel.scatter(kListF2good,EF2good[:,i],c=mF2good[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.',lw=0)
panel.set_xlim(-1.0,1.0)
panel.set_xlabel(r'$q$ [$k_L$]')
panel.set_title('(b)')

#make colorbar for F=2
cbar_axes = fig2.add_subplot(gs[0,8])
cbar = fig2.colorbar(d,cax=cbar_axes,ticks=np.array([-2,-1,0,1,2]))
cbar_axes.yaxis.set_ticks_position('left')

plt.savefig('loadF2.pdf', transparent=True)
