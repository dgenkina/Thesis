# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 17:43:49 2018

@author: dng5
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import lineFit2

rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['axes.titlesize'] = 8


rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']

#rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage[cm]{sfmath}']

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


figure=plt.figure()
figure.set_size_inches(3.5,3.0)
pan=figure.add_subplot(111)

xlist = np.linspace(-1.0,1.0,200)
ylist1=-np.cos(1.3*xlist*(2.0*np.pi))+1.0
ylist2=-np.cos((1.3*xlist+1.0/3.0)*(2.0*np.pi)) + 6.0
pan.plot(xlist,ylist1,'g-',linewidth=2)
pan.plot(xlist,ylist2,'b-',linewidth=2)
pan.set_ylim(-4.0,7.5)
pan.spines['top'].set_visible(False)
pan.spines['right'].set_visible(False)
pan.spines['bottom'].set_visible(False)
pan.spines['left'].set_visible(False)

pan.set_xlabel(r'Crystal momentum $q_x$')
pan.set_ylabel(r'Energy')
pan.set_xticks([])
pan.set_yticks([])

plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure3v3c.pdf',dpi=400)

figure=plt.figure()
figure.set_size_inches(3.5,6.0)
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.95, top=0.1, bottom = 0.05)
pan=figure.add_subplot(gs[0])
pan.set_title('3 sites \n 5 sites \n \\textbf{a} \n \\textbf{b} \n \\textbf{c} \n \\textbf{d} \n \\textbf{e} \n \\textbf{f} \n $\Phi/\Phi_0\\approx4/3$ \n $\Phi/\Phi_0=0$ \n $\\Phi/\Phi_0\\approx-4/3$ \n site $m$'
                '\n $m=-2$ \n $m=-1$ \n $m=0$ \n $m=1$ \n $m=2$ \n $j=-1$ \n $j=0$ \n $j=1$ \n $j=2$ \n'
                '$\phi=$ \n $0$ \n $2\pi/3$ \n $-2\pi/3$ \n $2k_R$ \n $2k_L$')
#                '\n Raman transition \n Reciprocal lattice period \n Fractional population'
#                '\n $F_x$ \n $F_{\parallel}$ \n Force \n Time \n $t$ \n $1$ \n $2$ \n $-1$ \n $-2$ \n $0$'
#                '\n $\mathbf{e_x}$ \n $\mathbf{e_s}$ \n $t_s$ \n $|t_x|e^{i\phi m}$'
#                '\n  $\\bar{m}=1$ \n  $\\bar{m}=0$ \n  $\\bar{m}=-1$'
#                '\n $q_s$ [$2k_L$] \n Energy \n $2k_L/Q$ \n $\Delta t$ \n $\Delta q_x$')
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/text.pdf',dpi=400)