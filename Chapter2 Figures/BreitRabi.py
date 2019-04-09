# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:11:32 2019

@author: dng5
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rcParams


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
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage[cm]{sfmath}',r'\usepackage{braket}']
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

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

gI = -0.000995141
gJ = 2.002331
h = 6.626068*1e-34
Ahfs = h*3.417*1e9
nuc = 3.0/2.0
ehfs = Ahfs*(nuc + 1.0/2.0)
mu = 9.274 * 1e-24
gImu = gI*mu/ehfs #in units of ehfs
gJmu = gJ*mu/ehfs # in units of ehfs


def x(B):
    return (gJmu - gImu)*B

def E1(Blist, mF): # in units of ehfs
    xlist = x(Blist)
    E = -1.0/(2.0*(2.0*nuc + 1.0)) + gImu*mF*Blist - (1.0/2.0)*np.sqrt(1.0 + (4.0*mF*xlist)/(2.0*nuc + 1.0) + xlist**2.0)
    return E
    
def E2(Blist, mF):
    xlist = x(Blist)
    E = -1.0/(2.0*(2.0*nuc + 1.0)) + gImu*mF*Blist + (1.0/2.0)*np.sqrt(1.0 + (4.0*mF*xlist)/(2.0*nuc + 1.0) + xlist**2.0)  
    return E
    
def E22(Blist, mF):
    E = 1.0*nuc/(2.0*nuc + 1.0) + np.float(np.sign(mF))*(gJmu + 2.0*nuc*gImu)*Blist/2.0
    return E

Blist = np.linspace(0.0,1.0,100)

figure = plt.figure(1)
figure.clear()
figure.set_size_inches(6.0,3.0)

gs = gridspec.GridSpec(1,2)
gs.update(left=0.1, right=0.95, top=0.85, bottom = 0.15, hspace=0.2,wspace=0.3)

pan=figure.add_subplot(gs[0])
for mF in np.arange(3)-1:
    pan.plot(Blist, E1(Blist,mF))
    pan.plot(Blist, E2(Blist,mF))
pan.plot(Blist, E22(Blist,2))
pan.plot(Blist, E22(Blist,-2))
pan.set_xlabel(r'$B$ [Tesla]')
pan.set_ylabel(r'Energy [$\Delta E_{\rm hfs}$]')
pan.set_title('(a)')
pan.text(0.0,0.8,r'$F=2$')
pan.text(0.0,-1.2,r'$F=1$')
pan.text(0.6,2.3,r'$m_J=+1/2$')
pan.text(0.6,-2.6,r'$m_J=-1/2$')

gI = 0.00017649
gJ = 2.0023193
Ahfs = -h*285.7308*1e6
nuc = 4.0
ehfs = Ahfs*(nuc + 1.0/2.0);
gImu = gI*mu/ehfs #in units of ehfs
gJmu = gJ*mu/ehfs # in units of ehfs


Blist = np.linspace(0.0,0.2,100)
pan=figure.add_subplot(gs[1])
for mF in np.arange(8)-3.5:
    pan.plot(Blist, -E1(Blist,mF)) #minus signs on the energies because ehfs is negative
    pan.plot(Blist, -E2(Blist,mF))
pan.plot(Blist, -E22(Blist,9.0/2.0))
pan.plot(Blist, -E22(Blist,-9.0/2.0))
pan.set_xlabel(r'$B$ [Tesla]')
pan.set_ylabel(r'Energy [$\Delta E_{\rm hfs}$]')
pan.set_title('(b)')
pan.text(0.0,1.1,r'$F=7/2$')
pan.text(0.0,-1.3,r'$F=9/2$')
pan.text(0.13,2.5,r'$m_J=+1/2$')
pan.text(0.13,-2.8,r'$m_J=-1/2$')

plt.savefig('BreitRabi.pdf',transparent=True)