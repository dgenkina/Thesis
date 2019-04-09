# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 17:03:47 2019

@author: swooty
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

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


xlist = np.linspace(0.01,3.0,200)

def bound(x):
    return -100000.0/x + 12000.0/(x**2.0)+30000.0

def free(x):
    return -100000.0/x + 30000.0/(x**2.0)

plt.close(1)
fig = plt.figure(1,figsize=(6.0,2.5))
gs = gridspec.GridSpec(1,2)
pan = fig.add_subplot(gs[0])
pan.plot(xlist,free(xlist),'g-')
pan.plot(xlist,bound(xlist),'r-')
pan.axhline(y=0.0,xmin=0.135,linestyle='--',color='g')
pan.axhline(y=-11000.0,xmin=0.075,xmax=0.74,linestyle='--',color='r')
pan.text(0.4,4000.0,'open channel',color='g')
pan.text(0.35,-25000.0,'closed channel',color='r')
pan.set_ylim(-200000.0,100000.0)
pan.set_xticks([])
pan.set_yticks([])
pan.set_ylabel('Energy')
pan.set_xlabel('Interatomic distance')
pan.set_title('(a)')

pan = fig.add_subplot(gs[1])
xlist = np.linspace(-4.0,4.0,200)
ylist = 0.1*xlist
pan.plot(xlist,np.zeros(xlist.size),'r--')
pan.plot(xlist,ylist,'g--')
pan.axvline(x=0.0,linestyle='-',color='k')
pan.text(0.1,-0.9,r'$B_0$')
pan.set_ylim(-1.0,1.0)
pan.set_xticks([])
pan.set_yticks([])
pan.set_xlabel('Magnetic field')
pan.set_title('(b)')

plt.savefig('Feshbach_plots.pdf',transparent=True)