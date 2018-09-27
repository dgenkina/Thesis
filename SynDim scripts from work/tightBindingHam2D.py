# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 17:21:35 2018

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=3.0
c=4.0/q#0.0#1064.0/790.0#
periodic=False
Flat=True

def TBham2d(tx,tm,p,q,kj,km):

    H=np.zeros((q,q), dtype=complex)
    for i in np.arange(q):
        H[i,i]=-2.0*tm*np.cos(km+2.0*np.pi*np.mod(i*p,q)/q)
        j = np.mod(i + p, q)
        H[i,j] += -tx*np.exp(1.0j*kj)
        H[j,i] += -tx*np.exp(-1.0j*kj)

    return H
    
def getEigenspectrum(tx,tm,p,q, plot=True):
    kjList=np.linspace(-np.pi/q,np.pi/q,100)
    kmList=np.linspace(-np.pi,np.pi,100)
    Egrid=np.zeros((kjList.size,kmList.size))
    for j,kj in enumerate(kjList):
        for m, km in enumerate(kmList):
            H=TBham2d(tx,tm,p,q,kj,km)
            E,V = sLA.eigh(H)
            sort=np.argsort(E)
            Egrid[j,m] = E[sort][0]
    if plot:        
        figure=plt.figure()
        pan=figure.add_subplot(111)
        pan.imshow(Egrid, cmap='Greys', extent=(kmList[0],kmList[-1],kjList[0],kjList[-1]))
    return Egrid,kjList,kmList
    
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8

rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']


rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}',r'\usepackage{amsmath}']


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

Egrid0,kj0,km0 = getEigenspectrum(0.5,0.5,0,1, plot=False)
EgridQ3P1,kjQ3P1,kmQ3P1 = getEigenspectrum(0.5,0.5,1,3, plot=False)
EgridQ7P3,kjQ7P3,kmQ7P3 = getEigenspectrum(0.5,0.5,3,7, plot=False)

figure = plt.figure()
figure.clear()
figure.set_size_inches(1.8,2.5)

gs = gridspec.GridSpec(3,1, height_ratios=[21, 7,3])

gs.update(left=0.3, right=0.95, top=0.95, bottom = 0.2, hspace=0.1)

s=km0.size
pan0=figure.add_subplot(gs[0])
pan0.imshow(Egrid0, cmap='Greys',extent=(km0[0]/(2.0*np.pi),km0[-1]/(2.0*np.pi),kj0[0]/(2.0*np.pi),kj0[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_xticklabels([])
pan0.set_yticks([-0.5,-0.25,0.0,0.25,0.5])
#pan0.set_ylabel(r'$q_s$ [$2k_L$]')

pan1=figure.add_subplot(gs[1])
pan1.imshow(EgridQ3P1, cmap='Greys',extent=(kmQ3P1[0]/(2.0*np.pi),kmQ3P1[-1]/(2.0*np.pi),kjQ3P1[0]/(2.0*np.pi),kjQ3P1[-1]/(2.0*np.pi)))
pan1.set_xticks([-0.5,0.0,0.5])
pan1.set_xticklabels([])
pan1.set_yticks([-0.08,0.08])

#pan1.set_ylabel(r'$q_s$ [$2k_L$]')

pan2=figure.add_subplot(gs[2])
pan2.imshow(EgridQ7P3, cmap='Greys',extent=(km0[0]/(2.0*np.pi),km0[-1]/(2.0*np.pi),kjQ7P3[0]/(2.0*np.pi),kjQ7P3[-1]/(2.0*np.pi)))
pan2.set_xticks([-0.5,0.0,0.5])
pan2.set_yticks([-0.07,0.07])
#pan2.set_ylabel(r'$q_s$ [$2k_L$]')
pan2.set_xlabel(r'$q_x$ [$2k_L$]')
figure.text(0.05,0.6,r'$q_s$ [$2k_L$]', rotation='vertical', size=8)
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure3v4a.pdf', dpi=500)

figure = plt.figure()
figure.clear()
figure.set_size_inches(1.75,2.5)
gs = gridspec.GridSpec(3,1)
gs.update(left=0.18, right=0.97, top=0.95, bottom = 0.2, hspace=0.3)

x1=np.linspace(-3,3,100)
y1=np.sin((x1-3.0/2.0)*2.0*np.pi/6.0)
pan0=figure.add_subplot(gs[0])
pan0.plot(x1,y1,'g-')
pan0.set_xlim(-3,7)
pan0.set_xticks([0])
pan0.set_yticks([])

x2=np.linspace(-3,7,100)
y2=np.sin((x2-3.0/2.0)*2.0*np.pi/6.0)
pan1=figure.add_subplot(gs[1])
pan1.plot(x2,y2,'g-')
pan1.set_xlim(-3,7)
pan1.set_xticks([0,6])
pan1.set_xticklabels(['0',r's(2$k_L$)'])
pan1.set_yticks([])
pan1.set_ylabel('Energy')

y3=np.sin((x2+1.0/2.0)*2.0*np.pi/6.0)
pan2=figure.add_subplot(gs[2])
pan2.plot(x2,y2,'g-')
pan2.plot(x2,y3,'b-')
pan2.set_xlim(-3,7)
pan2.set_xticks([-2,6])
pan2.set_xticklabels([r's(2$k_L$)+C(2$k_R$)',r's(2$k_L$)'])
pan2.set_xlabel(r'$q_x$')
pan2.set_yticks([])
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure3v4b.pdf', dpi=500)