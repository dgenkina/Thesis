# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 10:54:40 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams

datafile=np.load('fig5data.npz')
I0range=datafile['I0range']
v=datafile['v']

v100=np.loadtxt('Mathematica_100us.dat')
v50=np.loadtxt('Mathematica_50us.dat')
v10=np.loadtxt('Mathematica_10us.dat')

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

velplot = plt.figure(1,figsize=(3.2,3.0))
velplot.clear()
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15)
velpanel=velplot.add_subplot(gs[0])
velpanel.plot(v10[:,1],v10[:,0], 'b-', linewidth=2, label=r'$t = 10 \mu \rm{s}$')
velpanel.plot(I0range, [v[i,30]for i in range(I0range.size)], 'bo',  ms=3)
velpanel.plot(v50[:,1],v50[:,0], 'g-', linewidth=2, label=r'$t = 50 \mu \rm{s}$')
velpanel.plot(I0range, [v[i,150]for i in range(I0range.size)], 'go',  ms=3)
velpanel.plot(v100[:,1],v100[:,0], 'r-', linewidth=2, label=r'$t = 100 \mu \rm{s}$')
velpanel.plot(I0range, [v[i,300]for i in range(I0range.size)], 'ro',  ms=3)
velpanel.set_xlabel(r'$I_0/I_{\rm{sat}}$')
velpanel.set_ylabel('v[m/s]')

handles, labels = velpanel.get_legend_handles_labels()
velpanel.text(7.7,2.3,labels[0], rotation=0, size=10)
velpanel.text(7.7,8.3,labels[1], rotation=8, size=10)
velpanel.text(7.7,12.3,labels[2], rotation=12, size=10)

plt.savefig('figure5.pdf', transparent=True)
velplot.show()