# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 15:29:36 2015

@author: dng5
"""

import matplotlib.pyplot as plt
#import numpy as np



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np

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

figure = plt.figure(1,figsize=(3.2,3.0))
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15)
panel = figure.add_subplot(gs[0])

file1=np.load('OD0.1data.npz')
file2=np.load('OD1.0data.npz')
file3=np.load('OD10.0data.npz')


panel.plot(file1['trange']*1e6, 1.0/file1['signalToNoise'],'b-', linewidth=2.0)
panel.plot(file2['trange']*1e6, 1.0/file2['signalToNoise'],'g-', linewidth=2.0)
panel.plot(file3['trange']*1e6, 1.0/file3['signalToNoise'],'r-', linewidth=2.0)
panel.set_yscale('log')
panel.set_ylabel(r'SNR')
panel.set_xlim(0,80)
panel.set_ylim(0.2,4000)
panel.set_xlabel(r'time [$\mu \rm{s}$]')

panel.text(60,6,r'$\sigma_0 n$=0.1',  size=10)
panel.text(60,44,r'$\sigma_0 n$=1.0', size=10)
panel.text(60,1510,r'$\sigma_0 n$=10',  size=10)

plt.savefig('figure9a.pdf', transparent=True)