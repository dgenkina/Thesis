# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 12:19:15 2019

@author: swooty
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

chernOfSize = np.load('ChernOfSystemSizeOmega05.npz')
chernOfOmegaNF = np.load('ChernOfOmegaNotFlat.npz')
chernOfOmegaF = np.load('ChernOfOmegaFlat.npz')


plt.close(1)
fig = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,2)
gs.update(left=0.15, right=0.95, top=0.85, bottom = 0.15,wspace = 0.2,hspace=0.2)
pan=fig.add_subplot(gs[0,0])
pan.errorbar(2*chernOfSize['Slist']+1,chernOfSize['chernList'],yerr=chernOfSize['delChernList'],fmt='bo')
pan.set_xlabel(r'Number of sites along $\textbf{e}_s$')
pan.set_ylabel(r'Chern number')
pan.set_title(r'(a)')
pan.axhline(y=1.0,color='k',linestyle='--')

pan=fig.add_subplot(gs[0,1])
pan.errorbar(chernOfOmegaF['omegaList'],chernOfOmegaF['chernList1'],yerr=chernOfOmegaF['delChernList1'],fmt='go', label = '3 sites')
pan.errorbar(chernOfOmegaF['omegaList'][1:],chernOfOmegaF['chernList2'][1:],yerr=chernOfOmegaF['delChernList2'][1:],fmt='ro', label = '5 sites')
pan.set_xlabel(r'Coupling strength $\hbar\Omega$ [$E_L$]')
pan.set_title(r'(b)')
#pan.set_ylabel(r'Chern number')
pan.axhline(y=1.0,color='k',linestyle='--')
plt.legend()

plt.savefig('ChernDependenceFig.pdf',transparent=True)
plt.savefig('ChernDependenceFig.png',transparent=True, dpi=400)