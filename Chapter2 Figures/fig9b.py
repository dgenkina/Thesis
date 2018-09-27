# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 16:39:38 2015

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

file1=np.load('OD0.1OfInotData.npz')
file2=np.load('OD1.0OfInotData.npz')
file3=np.load('OD10.0OfInotData.npz')


panel.plot(file1['I0range'], 1.0/file1['signalToNoise'], 'b-',linewidth=2.0)
panel.plot(file2['I0range'], 1.0/file2['signalToNoise'],'g-', linewidth=2.0)
panel.plot(file3['I0range'], 1.0/file3['signalToNoise'], 'r-',linewidth=2.0)
panel.set_yscale('log')
panel.set_ylabel(r'SNR')
panel.set_xlabel(r'$I_0/I_{sat}$')
panel.set_ylim(4.0,2000)

panel.text(4,5.8,r'$\sigma_0 n$=0.1',  size=10)
panel.text(4,43,r'$\sigma_0 n$=1.0', size=10)
panel.text(4,1350,r'$\sigma_0 n$=10',  size=10)

plt.savefig('figure9b.pdf', transparent=True)