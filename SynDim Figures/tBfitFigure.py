# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 14:43:12 2018

@author: swooty
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib import rc, rcParams
import readIgor

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


#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,2)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.2,hspace=0.1)

F1file = np.load('TBfitToFullF1.npz')
panel=fig.add_subplot(gs[0])
lw=5
d=panel.plot(F1file['kList'],(F1file['Etb'][:,0]-F1file['offset2']+F1file['offset1']),'b-',lw=lw,label = 'Tight binding') #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F1file['kList'],F1file['Efull'][:,0],'r-',lw=lw, label = 'Full Hamiltonian') #c=mFull[:,i],vmin=-S,vmax=S,
d=panel.plot(F1file['kList'],(F1file['Etb'][:,1]-F1file['offset2']+F1file['offset1']),'b-',lw=lw) #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F1file['kList'],F1file['Efull'][:,1],'r-',lw=lw) #c=mFull[:,i],vmin=-S,vmax=S,
d=panel.plot(F1file['kList'],(F1file['Etb'][:,2]-F1file['offset2']+F1file['offset1']),'b-',lw=lw) #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F1file['kList'],F1file['Efull'][:,2],'r-',lw=lw)
panel.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')   
panel.set_ylabel(r'Energy [$E_L$]')
plt.legend(loc='upper right')
panel.set_title('(a)')
print F1file['ts']

F2file = np.load('TBfitToFullF2.npz')
panel=fig.add_subplot(gs[1])
lw=5
d=panel.plot(F2file['kList'],(F2file['Etb'][:,0]-F2file['offset2']+F2file['offset1']),'b-',lw=lw,label = 'Tight binding') #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F2file['kList'],F2file['Efull'][:,0],'r-',lw=lw, label = 'Full Hamiltonian') #c=mFull[:,i],vmin=-S,vmax=S,
d=panel.plot(F2file['kList'],(F2file['Etb'][:,1]-F2file['offset2']+F2file['offset1']),'b-',lw=lw) #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F2file['kList'],F2file['Efull'][:,1],'r-',lw=lw) #c=mFull[:,i],vmin=-S,vmax=S,
d=panel.plot(F2file['kList'],(F2file['Etb'][:,2]-F2file['offset2']+F2file['offset1']),'b-',lw=lw) #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(F2file['kList'],F2file['Efull'][:,2],'r-',lw=lw)
panel.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')
panel.set_title('(b)')   
#panel.set_ylabel(r'Energy [$t_x$]')
plt.legend(loc='upper right')

print F2file['ts']

plt.savefig('tbFit.pdf',transparent=True)