# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 21:54:59 2019

@author: swooty
"""
import numpy as np
from scipy import linalg as sLA
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec


rcParams['axes.labelsize'] = 10
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8

rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}', r'\usepackage{braket}']

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

hofFile = np.load('HofstadtersButt.npz')
flux = hofFile['fluxList']
E = hofFile['ElistAll']/0.1

E = E[flux!=1]
flux = flux[flux!=1]

plt.close(1)
fig = plt.figure(1,figsize=(6.2,5.0))
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15,hspace = 0.3,wspace=4.5)


pan = fig.add_subplot(gs[0])
#pan.plot(hofFile['fluxList'],hofFile['ElistAll']/0.1,'k.', markersize=1)
pan.hist2d(flux,E,bins=[620,400],normed=True,vmax=2, cmap='gray_r')
pan.set_xlabel(r'Flux $\Phi/\Phi_0$')
pan.set_xlim(0.0,1.0)
pan.set_ylabel(r'Energy E/$t_x$')

#plt.savefig('HofstadterFig.pdf',transparent=True)