# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 17:05:53 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams


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

datafile = np.load('fig7data.npz')
I0range=datafile['I0range']
nu1=datafile['nu1']
nu2=datafile['nu2']
nuS=datafile['nuS']
nuM=datafile['nuM']
tind = 299

nuplot = plt.figure()
nuplot.clear()
nuplot.set_size_inches(3.2,3.0)
gs = gridspec.GridSpec(2,1)
gs.update(left=0.2, right=0.95, top=0.95, bottom = 0.15)
nupanel=nuplot.add_subplot(gs[0])
nupanel.plot(I0range, [nu1[i,tind] for i in range(I0range.size)], 'k-', linewidth=1, label=r'$OD_1$')
nupanel.plot(I0range, [nu2[i,tind] for i in range(I0range.size)], 'b-', linewidth=1, label=r'$OD_2$')
nupanel.plot(I0range, [nuS[i,tind] for i in range(I0range.size)], 'g-', linewidth=1, label=r'$OD_{corr1}$')
nupanel.plot(I0range, [nuM[i,tind] for i in range(I0range.size)], 'go', ms=2, label=r'$OD_{corr2}$')
nupanel.axhline(y=0,linestyle='--',color='k')
nupanel.set_ylabel('Optical Depth')
nupanel.set_xticks([])

nupanel.text(3.5,1.1,r'$OD^{(1)}$', size=10)
nupanel.text(3.5,-2.5,r'$OD^{(2)}$', size=10)
nupanel.text(3.5,0.27,r'$OD^{\rm{sim}}$',  size=10)

diffpanel=nuplot.add_subplot(gs[1])
diffpanel.plot(I0range,[(nuS[i,tind]-nuM[i,tind])/nuM[i,tind] for i in range(I0range.size)],'ko')
diffpanel.set_ylim([0.0,0.005])
diffpanel.set_xlabel(r'$I_0/I_{\rm sat}$')
diffpanel.set_ylabel('Fractional difference')
nuplot.show()
plt.savefig('figure7.pdf', transparent=True)