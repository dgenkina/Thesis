# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 15:50:33 2019

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
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{braket}']

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


hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg


A = 0.7e6 #MHz per Gauss, approximate for F=1 of Rb
omega_rf = 1.0e6 #rf angular frequency

def delta(B):
    return 2.0*np.pi*A*B - omega_rf

def E(B, Omega):
    E = np.sqrt(delta(B)**2.0 + Omega**2.0)/2
    return (-E,E)

Blist = np.arange(0.0,0.455,0.0001)
EoffP = np.zeros(Blist.size)
EoffM = np.zeros(Blist.size)
EonP = np.zeros(Blist.size)
EonM = np.zeros(Blist.size)
Omega = 2.0e5

for ind,B in enumerate(Blist):
    EoffM[ind] = E(B,0.0)[0]
    EoffP[ind] = E(B,0.0)[1]
    EonM[ind] = E(B, Omega)[0]
    EonP[ind] = E(B,Omega)[1]

#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.0,3.0))
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.5,hspace=1.0)
pan = fig.add_subplot(gs[0])

pan.plot(Blist,EoffM*1e-5, 'k-')
pan.plot(Blist,EoffP*1e-5, 'k-')
pan.plot(Blist,EonM*1e-5, 'b-')
pan.plot(Blist,EonP*1e-5, 'b-')
pan.text(0,-4.1,r'$\ket{g}$')
pan.text(0,4.0, r'$\ket{e}$')
pan.set_ylabel(r'Energy [arb. u.]')
pan.set_xlabel(r'B [Gauss]')
pan.axvline(x=0.0275, color='k',linestyle= '--')
pan.axvline(x=Blist[Blist.size/2],color='k',linestyle= '--')
pan.axvline(x=0.4275, color='k',linestyle= '--')
pan.text(0.03, -5.4, r'$B_{start}$')
pan.text(0.23, -5.4, r'$B_{res}$')
pan.text(0.385, -5.4, r'$B_{stop}$')
plt.savefig('ARPfigure.pdf',transparent=True)