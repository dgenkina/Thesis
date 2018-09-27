# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 10:32:20 2014

@author: dng5
"""

import matplotlib.pyplot as plt
#import numpy as np
import ImagingModel3Fig1



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np

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

# Pull data from hdf5 file
#data1 = h5scripting.get_all_data(h5name, groupname1)
#data2 = h5scripting.get_all_data(h5name, groupname2)

#boxpoints = np.array([
#                        [2,2.25],
#                        [2,4.5],
#                        [6,4.5],
#                        [6,2.25],
#                        [2,2.25]
#                     ])

# Make the actual plot

tfinal = 0.00005
steps = 300
trange = np.linspace(0,tfinal,steps)
outputTuple = ImagingModel3Fig1.Image(1.5,tfinal,steps,1.6)
superAtomNumber = len(outputTuple[5])
dz = outputTuple[6]

ods = plt.figure(1,figsize=(3.2,3.0))
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15)
#odplot, ax1 = plt.subplots()
ax1=ods.add_subplot(gs[0])
#tind = 0
#line0 = ax1.plot(np.arange(superAtomNumber)*dz*1e6, [outputTuple[5][i][tind] for i in range(superAtomNumber)], 'y-', linewidth=2, label="t="+str(int(trange[tind]*1e6))+" $\mu s$")
tind = 60
line1 = ax1.plot(np.arange(superAtomNumber)*dz*1e6, [outputTuple[5][i][tind] for i in range(superAtomNumber)], 'r-', linewidth=2, label="t="+str(int(trange[tind]*1e6))+" $\mu s$")
tind = 150
line2 = ax1.plot(np.arange(superAtomNumber)*dz*1e6, [outputTuple[5][i][tind] for i in range(superAtomNumber)], 'b-', linewidth=2, label="t="+str(int(trange[tind]*1e6))+" $\mu s$")
tind = 299
line3 = ax1.plot(np.arange(superAtomNumber)*dz*1e6, [outputTuple[5][i][tind] for i in range(superAtomNumber)], 'g-', linewidth=2, label="t="+str(int(trange[tind]*1e6))+" $\mu s$")
ax1.set_xlabel("z [$\mu m$]")
ax1.set_ylabel("Velocity [$m/s$]")
#ax1.legend(loc=3)

ax2=ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[0],ax1.get_ylim()[1]*k*2.0/Gamma)
ax2.set_ylabel("$2\delta /\Gamma$")

handles, labels = ax1.get_legend_handles_labels()
ax1.text(0.375,1.25,labels[0], rotation=-7, size=10)
ax1.text(0.375,2.65,labels[1], rotation=-7, size=10)
ax1.text(0.375,4.33,labels[2], rotation=-8, size=10)

ods.show()
plt.savefig('figure1.pdf', transparent=True)


#fig=mpl.pyplot.figure(1,figsize=(6.0,2.75))
#gs = gridspec.GridSpec(1, 2)
#gs.update(left=0.1, right=0.95, top=0.92, bottom = 0.15, hspace=0.3, wspace = 0.25)
#
#ax=fig.add_subplot(gs[0])
#
#ax.plot(data1["Omegas"], data1["Energies"],"-", zorder=2)
#ax.plot(boxpoints[:,0],boxpoints[:,1], "-", color="gray", linewidth=2, zorder=1)
#
#ax.set_xlabel('Raman coupling in units of $E_L$')
#ax.set_ylabel('Energy in units of $E_L$')
#ax.set_ylim(-5,10)
#
#ax.set_title(r'(a) Exact energies with harmonic confinment', fontsize=8, loc = 'left')
#
#ax=fig.add_subplot(gs[1])
#
#ax.plot(data2["Omegas"], data2["Energies"],"-")
#
#ax.set_xlabel('Raman coupling in units of $E_L$')
#ax.set_ylabel('Energy in units of $E_L$')
#ax.set_ylim(2.25,4.5)
#ax.set_xlim(2,6)
#
#ax.set_title(r'(b) Zoomed view', fontsize=8, loc = 'left')
#
#
#if pdfName is not None:
#    mpl.pyplot.savefig(pdfName, transparent=True)
#
#mpl.pyplot.show()
#mpl.pyplot.clf()


