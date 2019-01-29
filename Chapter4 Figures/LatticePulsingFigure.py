# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 15:48:35 2018

@author: swooty
"""

import numpy as np
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib import rc, rcParams
import readIgor

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
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{upgreek}']

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


#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,4.0))
gs = gridspec.GridSpec(2,3)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.5,hspace=0.5)
 
filename = 'C:/Users/swooty/Documents/Thesis Data/2017Aug20 Lattice calibrations/PIXIS_20Aug2017_0066.ibw' 

od = readIgor.processIBW(filename)['rotOD'][425:615,445:520]
pan = fig.add_subplot(gs[:,0])
pan.imshow(od,vmin=-0.1, vmax=0.3)
pan.set_xticks([])
pan.set_yticks([32,97,162])
pan.set_yticklabels([r'$2k_{\rm L}$',r'$0$',r'$-2k_{\rm L}$'])
pan.set_title('(a)')
pan.set_ylabel(r'$k_x$')

# Pulsing panel for lower coupling 
lowU = np.load('LatticePulsingLowU.npz')

panel = fig.add_subplot(gs[0,1:])
panel.plot(lowU['time']*1e6,lowU['fractions'][:,0],'bo', label=r'$k_x=+2k_{\rm R}$')
panel.plot(lowU['time']*1e6,lowU['fractions'][:,1],'go', label=r'$k_x=0$')
panel.plot(lowU['time']*1e6,lowU['fractions'][:,2],'ro', label=r'$k_x=-2k_{\rm R}$')
panel.plot(lowU['tForFit']*1e6,lowU['pop0'],'b-')
panel.plot(lowU['tForFit']*1e6,lowU['pop1'],'g-')
panel.plot(lowU['tForFit']*1e6,lowU['pop2'],'r-')

panel.set_ylabel('Fractional population')
panel.set_title('(b)')


# Pulsing panel for higher coupling 
highU = np.load('LatticePulsingHighU.npz')

panel = fig.add_subplot(gs[1,1:])
panel.plot(highU['time']*1e6,highU['fractions'][:,0],'bo', label=r'$k_x=+2k_{\rm R}$')
panel.plot(highU['time']*1e6,highU['fractions'][:,1],'go', label=r'$k_x=0$')
panel.plot(highU['time']*1e6,highU['fractions'][:,2],'ro', label=r'$k_x=-2k_{\rm R}$')
panel.plot(highU['tForFit']*1e6,highU['pop0'],'b-')
panel.plot(highU['tForFit']*1e6,highU['pop1'],'g-')
panel.plot(highU['tForFit']*1e6,highU['pop2'],'r-')

panel.set_xlabel(r'Lattice pulse time [$\upmu$s]')
panel.set_ylabel('Fractional population')
panel.set_title('(c)')
plt.legend()


fig.show()

plt.savefig('LatticePulsing.pdf', transparent=True)