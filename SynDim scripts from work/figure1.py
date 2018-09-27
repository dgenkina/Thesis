# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:20:57 2017

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

datafile1=np.load('SynDimBandStructure_F10_n15_Chern1Flat.npz')
datafile3=np.load('SynDimBandStructure_F1_n7_Chern1.npz')
datafile2=np.load('SynDimBandStructure_F2_n7_Chern1.npz')
datafile4=np.load('SynDimBandStructure_F2_n7_Chern0.npz')
datafile5=np.load('SynDimBandStructure_F1_n7_Chern1_ExtendedFlat.npz')
lw=4
width=1

fig = plt.figure()
fig.clear()
fig.set_size_inches(3.7,5.0)
gs = gridspec.GridSpec(3,2)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15)
gs.update(hspace=0.5,wspace=0.5)

#ax0=fig.add_subplot(gs[0])
#ax0.imshow(np.load('22Aug2017_0006_processedImage.npz')['odFiltered'])


ax1=fig.add_subplot(gs[2])
S1=10
for i in range(2*S1+1):   
    d=ax1.scatter(datafile1['kList'],datafile1['E'][:,i],c=datafile1['m'][:,i],vmin=-S1,vmax=S1, marker='_')
plt.colorbar(d)   
ax1.set_ylabel(r'Energy [$E_L$]') 
ax1.set_xlabel(r'$q_x$ [$k_L$]')
ax1.xaxis.set_ticks([-1.0,0.0,1.0])
ax1.set_title('(b)',loc='left')

#
#ax2=fig.add_subplot(gs[3])      
#for i in range(2*S1+1):    
#    ax2.scatter(datafile1['kList'],[i-S1 for j in range(datafile1['kList'].size)],c=datafile1['popAvg'][:,i],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S1,marker='_',linewidths=lw)
#ax2.set_ylim(-S1-1,S1+1)
#ax2.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
#ax2.set_ylabel(r'Lattice site in $\hat{y}$')

ax3=fig.add_subplot(gs[1])
S=2
for i in range(2*S+1):   
    d=ax3.scatter(datafile2['kList'],datafile2['E'][:,i],c=datafile2['m'][:,i],vmin=-S,vmax=S, marker='_')
plt.colorbar(d)  
#ax3.set_ylabel(r'Energy [$E_L$]') 
ax3.xaxis.set_ticks([])
ax3.set_title('(c)',loc='left')


#ax4=fig.add_subplot(gs[4])      
#for i in range(2*S+1):    
#    ax4.scatter(datafile2['kList'],[i-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,i],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
#ax4.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
#ax4.set_ylim(-S1-1,S1+1)
#ax4.yaxis.set_ticks([])

ax4=fig.add_subplot(gs[0])
S=1
for i in range(2*S+1):   
    d=ax4.scatter(datafile3['kList'],datafile3['E'][:,i],c=datafile3['m'][:,i],vmin=-S,vmax=S, marker='_')
plt.colorbar(d)    
ax4.set_ylabel(r'Energy [$E_L$]') 
ax4.xaxis.set_ticks([])
ax4.set_title('(a)',loc='left')


ax5=fig.add_subplot(gs[3])
S=2
for i in range(2*S+1):   
    d=ax5.scatter(datafile4['kList'],datafile4['E'][:,i],c=datafile4['m'][:,i],vmin=-S,vmax=S, marker='_')
plt.colorbar(d)    
#ax5.set_ylabel(r'Energy [$E_L$]')     
ax5.set_xlabel(r'$q_x$ [$k_L$]')
ax5.set_title('(d)',loc='left')
ax5.xaxis.set_ticks([-1.0,0.0,1.0])



ax6=fig.add_subplot(gs[4:6])
S=1
for i in range(2*S+1):   
    d=ax6.scatter(datafile5['kList'],datafile5['E'][:,i],c=datafile5['m'][:,i],vmin=-S,vmax=S, marker='_')
#colorbar(d)    
ax6.set_ylabel(r'Energy [$E_L$]')     
ax6.set_xlabel(r'$q_x$ [$k_L$]')
ax6.yaxis.set_ticks([-0.6,-0.4,-0.2])
ax6.set_title('(e)',loc='left')

#ax6=fig.add_subplot(gs[5])      
#for i in range(2*S+1):    
#    ax6.scatter(datafile3['kList'],[i-S for j in range(datafile3['kList'].size)],c=datafile3['pops'][:,0,i],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
#ax6.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
#ax6.set_ylim(-S1-1,S1+1)
#ax6.yaxis.set_ticks([])
#plt.tight_layout(pad=1.08,h_pad=0.05,w_pad=None)
fig.show()
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/theoryPlots.pdf', transparent=True)