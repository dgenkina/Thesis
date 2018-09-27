# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 18:31:50 2018

@author: swooty
"""

import numpy as np
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle

import fitLattice

from matplotlib import rc, rcParams


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

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

def plotLatticeBandStructure(U,n, BZs = 1, plot = True):
    kList = np.linspace(-BZs,BZs,600)
    Nlat=2*n+1
    Energies = np.zeros((kList.size,Nlat))
    for ind, k in enumerate(kList):
        H = fitLattice.LatHam(U,n,k)
        Energy, V = sLA.eig(H)
        sort = np.argsort(Energy)
        Energies[ind] = Energy[sort]
        if ind == kList.size/2:
            eigV = V[:,sort]
        
    print np.max(Energies[:,0])-np.min(Energies[:,0])
        
    if plot:    
        figure = plt.figure()
        pan = figure.add_subplot(111)
        pan.plot(kList, Energies[:,:3], 'b.')
        pan.set_xlabel(r'$q/k_L$')
        pan.set_ylabel('Energy')
    
    return kList, Energies, eigV

kList0, Energies0, eigV = plotLatticeBandStructure(0,7, BZs = 3, plot = False)
kList4, Energies4, eigV = plotLatticeBandStructure(4.0,7, BZs = 3, plot = False)

#figure = plt.figure(1,figsize=(6.2,4.0))
#gs = gridspec.GridSpec(1,1)
#gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15)
#pan1 = figure.add_subplot(gs[0])
#pan1.plot(kList0, Energies0[:,:3], 'b-.')
#pan1.plot(kList4,Energies4[:,:3], 'b-')
#pan1.set_xlabel(r'$q/k_L$')
#pan1.set_ylabel(r'Energy [$E_L$]')
#pan1.axvline(x=1.0,c='k')
#pan1.axvline(x=-1.0,c='k')
#
#plt.savefig('figure1.pdf', transparent=True)
#figure.show()
#    
#
#pops = fitLattice.propagateLatHamFit(0.00012*Erecoil/hbar, 8.0,k=0,psi0='norm' )
#
#
kList0, Energies0, eigV0 = plotLatticeBandStructure(0,7, BZs = 1, plot = False)
kList4, Energies4, eigV4 = plotLatticeBandStructure(8.0,7, BZs = 1, plot = False)
#
#psi0 = np.zeros(Energies4.shape[1],dtype=complex)
#psi0[psi0.size/2] = 1.0
#
#Vinv4 = sLA.inv(eigV4)
#psi4 = np.dot(Vinv4,psi0)
#proj4 = np.absolute(psi4)**2.0
#
#Vinv0 = sLA.inv(eigV0)
#psi0 = np.dot(Vinv0,psi0)
#proj0 = np.absolute(psi0)**2.0
#
#scale = 20
#
#fig2 = plt.figure(2,figsize=(6.2,4.0))
#gs = gridspec.GridSpec(4,3)
#gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,hspace = 1.0)
#
#pan = fig2.add_subplot(gs[0,:])
#xlist = np.linspace(0,10,1000)
#ylist = 8.0*np.ones(xlist.size)
#ylist[xlist<2.0] = 0.0
#ylist[xlist>8.0] = 0.0
#pan.plot(xlist,ylist,'k-')
#pan.set_xticks([])
#pan.set_xlabel(r'time')
#pan.set_yticks([0,8])
#pan.set_ylabel(r'$V_0$ [$E_L$]')
#pan.set_title('(a)')
#
#
#
#pan = fig2.add_subplot(gs[1:,1])
#pan.plot(kList4,Energies4[:,:3], 'b-') 
#pan.plot(0,Energies4[kList4.size/2,0], marker = 'o',c='b', markersize = np.sqrt(proj4[0])*scale)
#pan.plot(0,Energies4[kList4.size/2,1], marker = 'o',c='b', markersize = np.sqrt(proj4[1])*scale)
#pan.plot(0,Energies4[kList4.size/2,2], marker = 'o',c='b', markersize = np.sqrt(proj4[2])*scale)
#pan.set_yticks([])
#pan.set_xlabel(r'$q/k_L$')
#pan.set_title('(c)')
#
#pan = fig2.add_subplot(gs[1:,0])
#pan.plot(kList0,Energies0[:,:3], 'b-') 
#pan.plot(0,Energies0[kList0.size/2,0], marker = 'o',c='b', markersize = np.sqrt(proj0[0])*scale)
#pan.plot(0,Energies0[kList0.size/2,1], marker = 'o',c='b', markersize = np.sqrt(proj0[1])*scale)
#pan.plot(0,Energies0[kList0.size/2,2], marker = 'o',c='b', markersize = np.sqrt(proj0[2])*scale)
#pan.set_xlabel(r'$q/k_L$')
#pan.set_ylabel(r'Energy [$E_L$]')
#pan.set_title('(b)')
#
#pan = fig2.add_subplot(gs[1:,2])
#pan.plot(kList0,Energies0[:,:3], 'b-') 
#pan.plot(0,Energies0[kList0.size/2,0], marker = 'o',c='b', markersize = np.sqrt(pops[1])*scale)
#pan.plot(0,Energies0[kList0.size/2,1], marker = 'o',c='r', markersize = np.sqrt(pops[0])*scale*np.sqrt(2),fillstyle='left',markeredgewidth = 0)
#pan.plot(0,Energies0[kList0.size/2,2], marker = 'o',c='g', markersize = np.sqrt(pops[2])*scale*np.sqrt(2),fillstyle = 'right', markeredgewidth = 0)
#pan.set_yticks([])
#pan.set_xlabel(r'$q/k_L$')
#pan.set_title('(d)')
#plt.savefig('figure2.pdf', transparent=True)    
#
#fig2.show()
#
#
#


#Adiabatic loading figure
pops = np.abs(eigV4[:,0])**2.0

fig2 = plt.figure(3,figsize=(6.2,4.0))
gs = gridspec.GridSpec(4,3)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,hspace = 1.0)

pan = fig2.add_subplot(gs[0,:])
xlist = np.linspace(0,10,1000)
ylist = (xlist - 2.0)*8.0/6.0
ylist[xlist<2.0] = 0.0
ylist[xlist>8.0] = 0.0
pan.plot(xlist,ylist,'k-')
pan.set_xticks([])
pan.set_xlabel(r'time')
pan.set_yticks([0,8])
pan.set_ylabel(r'$V_0$ [$E_L$]')
pan.set_title('(a)')



pan = fig2.add_subplot(gs[1:,1])
pan.plot(kList4,Energies4[:,:3], 'b-') 
pan.plot(0,Energies4[kList4.size/2,0], marker = 'o',c='b', markersize = np.sqrt(proj0[0])*scale)
pan.plot(0,Energies4[kList4.size/2,1], marker = 'o',c='b', markersize = np.sqrt(proj0[1])*scale)
pan.plot(0,Energies4[kList4.size/2,2], marker = 'o',c='b', markersize = np.sqrt(proj0[2])*scale)
pan.set_yticks([])
pan.set_xlabel(r'$q/k_L$')
pan.set_title('(c)')

pan = fig2.add_subplot(gs[1:,0])
pan.plot(kList0,Energies0[:,:3], 'b-') 
pan.plot(0,Energies0[kList0.size/2,0], marker = 'o',c='b', markersize = np.sqrt(proj0[0])*scale)
pan.plot(0,Energies0[kList0.size/2,1], marker = 'o',c='b', markersize = np.sqrt(proj0[1])*scale)
pan.plot(0,Energies0[kList0.size/2,2], marker = 'o',c='b', markersize = np.sqrt(proj0[2])*scale)
pan.set_xlabel(r'$q/k_L$')
pan.set_ylabel(r'Energy [$E_L$]')
pan.set_title('(b)')

pan = fig2.add_subplot(gs[1:,2])
pan.plot(kList0,Energies0[:,:3], 'b-') 
pan.plot(0,Energies0[kList0.size/2,0], marker = 'o',c='b', markersize = np.sqrt(pops[Energies4.shape[1]/2])*scale)
pan.plot(0,Energies0[kList0.size/2,1], marker = 'o',c='r', markersize = np.sqrt(pops[Energies4.shape[1]/2-1])*scale*np.sqrt(2),fillstyle='left',markeredgewidth = 0)
pan.plot(0,Energies0[kList0.size/2,2], marker = 'o',c='g', markersize = np.sqrt(pops[Energies4.shape[1]/2+1])*scale*np.sqrt(2),fillstyle = 'right', markeredgewidth = 0)
pan.set_yticks([])
pan.set_xlabel(r'$q/k_L$')
pan.set_title('(d)')
plt.savefig('figure3.pdf', transparent=True)    

fig2.show()
        
 