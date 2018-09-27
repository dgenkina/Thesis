# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 14:58:30 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import ImagingModel3Fig4

tfinal = 0.00008
steps = 300
trange = np.linspace(0,tfinal,steps)
outputTuple = ImagingModel3Fig4.Image(1.2,tfinal,steps,1.6)
superAtomNumber = len(outputTuple[5])

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

distances = plt.figure(1,figsize=(3.2,3.0))
distances.clear()
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15)
dist=distances.add_subplot(gs[0])
dist.plot(1e6*trange, [1e6*outputTuple[5][0][i] for i in range(steps)], 'g-', linewidth=2, label="front")
dist.plot(1e6*trange, [1e6*outputTuple[5][(superAtomNumber-1)/2][i] for i in range(steps)], 'r-', linewidth=2, label="middle")
dist.plot(1e6*trange, [1e6*outputTuple[5][superAtomNumber-1][i] for i in range(steps)], 'b-', linewidth=2, label="back")
dist.set_ylabel(r'z position [$\mu \rm{m}$]')
dist.set_xlabel(r'time [$\mu \rm{s}$]')
dist.text(0,8,"front", rotation=8, size=10)
dist.text(0,35,"middle", rotation=8, size=10)
dist.text(0,56,"back", rotation=5, size=10)

plt.savefig('figure4.pdf', transparent=True)
distances.show()