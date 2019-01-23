# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 18:25:32 2014

@author: dng5
"""
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import ImagingModelFig2


#get the data
tfinal = 0.00005
steps = 300
trange = np.linspace(0,tfinal,steps)
outputTuple = ImagingModelFig2.Image(1.5,tfinal,steps,1.6)


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


tfinal = 0.0001
steps = 300
trange = np.linspace(0,tfinal,steps)
outputTuple = ImagingModelFig2.Image(0.8,tfinal,steps,1.6)

ods = plt.figure()
ods.clear()
ods.set_size_inches(3.2,3.0)
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15)
odpanel = ods.add_subplot(gs[0])
line0 = odpanel.plot(trange*1.0e6, [outputTuple[0] for i in trange], 'b-', linewidth=2, label="$\sigma_0 n$")
line1 = odpanel.plot(trange*1.0e6, outputTuple[2], 'g-', linewidth=2, label="$OD_1$")
line2 = odpanel.plot(trange*1.0e6, outputTuple[3], 'r-', linewidth=2, label="$OD_2$")
odpanel.set_xlabel("Imaging time [us]")
odpanel.set_ylabel("Optical Depth")
odpanel.axvline(x = 25.0,linestyle='--',color = 'k')

handles, labels = odpanel.get_legend_handles_labels()
odpanel.text(75,1.65,r'$\sigma_0 n$', rotation=0, size=10)
odpanel.text(75,1.0,'Eq. 2.27', rotation=-10, size=10)
odpanel.text(75,3.55,'Eq. 4.6', rotation=50, size=10)
odpanel.text(19,3.0, 'recoil time', rotation=90)

plt.savefig('figure2.pdf', transparent=True)
ods.show()