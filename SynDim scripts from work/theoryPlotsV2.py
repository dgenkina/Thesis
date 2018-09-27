# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:21:15 2017

@author: dng5
"""

import matplotlib.pyplot as plt
#import numpy as np



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np

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

ColorMap = "afmhot"


datafile1=np.load('SynDimBandStructure_F1_n7_Chern1_ExtendedFlat.npz')
datafile2=np.load('ChernOfOmegaNotFlat.npz')
datafile3=np.load('ChernOfSystemSize.npz')
lw=4
width=1

fig = plt.figure()
fig.clear()
fig.set_size_inches(3.7,5.0)
gs = gridspec.GridSpec(3,1)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15)
gs.update(hspace=0.5,wspace=0.5)


ax6=fig.add_subplot(gs[0])
S=1
for i in range(2*S+1):   
    d=ax6.scatter(datafile1['kList'],datafile1['E'][:,i],c=datafile1['m'][:,i],vmin=-S,vmax=S, marker='_')
#colorbar(d)    
ax6.set_ylabel(r'Energy [$E_L$]')     
ax6.set_xlabel(r'$q_x$ [$k_L$]')
ax6.yaxis.set_ticks([-0.6,-0.4,-0.2])
ax6.set_title('(a)',loc='left')



ax1=fig.add_subplot(gs[1])
S1=10
for i in range(2*S1+1):   
    d=ax1.plot(datafile2['omegaList'],datafile2['chernList1'],'b-')
    d=ax1.plot(datafile2['omegaList'][10:],datafile2['chernList2'][10:],'r-')
chernlistTot=np.append(datafile2['chernList1'], datafile2['chernList2'][10:])
ax1.set_xlabel(r'$\Omega$ [$E_L$]') 
ax1.set_ylabel(r'Chern number')
ax1.yaxis.set_ticks([np.round(np.min(chernlistTot),decimals=2),np.round(np.max(chernlistTot),decimals=2)])
ax1.set_title('(b)',loc='left')

ax2=fig.add_subplot(gs[2])
S1=10
for i in range(2*S1+1):   
    d=ax2.plot(datafile3['Slist']*2.0+1.0,datafile3['chernList'],'bo')
ax2.yaxis.set_ticks([np.round(np.min(datafile3['chernList']),decimals=2),np.round(np.max(datafile3['chernList']),decimals=2)])
ax2.set_xlabel(r'Strip width [sites]') 
ax2.set_ylabel(r'Chern number')

ax2.set_title('(c)',loc='left')

fig.show()
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/theoryPlotsV2.pdf', transparent=True)