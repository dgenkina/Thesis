# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 13:47:45 2014

@author: dng5
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np
from matplotlib.ticker import MultipleLocator, NullLocator, LogLocator

rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 10
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

datafile = np.load('fig3data.npz')
I0range = datafile['I0range']
nu=datafile['nu']
nu1=datafile['nu1']
trange=datafile['trange']
timeIndexes = np.array([18,150,299])   

odfigure = plt.figure()
odfigure.clear()
odfigure.set_size_inches(6,3)
gs = gridspec.GridSpec(1,3)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15, hspace=0.1, wspace=0.25)

# Look up add_subplots

odpanel={}
for t in range(len(timeIndexes)):
    index = timeIndexes[t]
    odpanel[t] = odfigure.add_subplot(gs[t])
    odpanel[t].plot(I0range, nu[:,index], 'bo', label=r'$\nu_{corr}$')
    odpanel[t].plot(I0range, nu1[:,index], 'g-',linewidth=2, label=r"$\nu_1$")
    odpanel[t].set_xlabel(r"$I_0/I_{\rm sat}$")
    odpanel[t].set_title(r"$t = {0:.0f}\ \mu{{\rm s}}$".format(trange[index]*1e6), fontsize=10)
    odpanel[t].set_xscale('log')
    odpanel[t].set_xlim(1e-2,1e4)
    odpanel[t].set_xticks([1e-2,1e1,1e4])
    if t>0:
        odpanel[t].yaxis.set_major_locator(NullLocator())

odpanel[0].set_ylabel(r"Optical Depth")




#odpanel[0].legend()
odfigure.show()
plt.savefig('figure3.pdf', transparent=True)