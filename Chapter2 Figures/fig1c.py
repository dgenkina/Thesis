# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 17:44:21 2015

@author: dng5
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import ImagingModel3

#generate the data
outputTuple = ImagingModel3.Image(1.2, 0.000001,2,0.00001)
zrange = outputTuple[8]
gaussian = outputTuple[5]
superAtomNumber = outputTuple[7]
atomZ = [zrange[0]*10e6]
[atomZ.append(outputTuple[6][i][0]*10e6) for i in range(superAtomNumber)]
atomZ.append(10e6*zrange[zrange.size-1])
probeProfile = [outputTuple[9][0,0]]
[probeProfile.append(outputTuple[9][0,i+1]) for i in range(superAtomNumber)]
probeProfile.append(outputTuple[9][0,superAtomNumber])

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


fig = plt.figure()
fig.clear
fig.set_size_inches(2.4,2.6)
gs = gridspec.GridSpec(2,1)
gs.update(left=0.2, right=0.95, top=0.85, bottom = 0.15)

#odplot, ax1 = plt.subplots()
panel1=fig.add_subplot(gs[0])
#panel1=figure1.add_subplot(2,1,1)
panel1.plot(zrange*10e6, gaussian, linewidth=2.0)
panel1.xaxis.set_ticks([])
panel1.yaxis.set_ticks(np.linspace(0, 3.4e7,3))
panel1.set_ylabel(r'n [atoms/$\rm{m}^2$]')
panel2=fig.add_subplot(gs[1])
panel2.plot(atomZ, probeProfile, 'k-', linewidth=2.0)
panel2.set_xlabel(r'z [$\mu \rm{m}$]')
panel2.set_ylabel(r'$I_0/I_{\rm{sat}}$')
panel2.yaxis.set_ticks(np.round(np.linspace(panel2.get_ylim()[0], panel2.get_ylim()[1],3),decimals=1))
panel2.set_ylim(0.4,1.3)

fig.savefig('Picture1c.pdf', transparent=True)
fig.show()