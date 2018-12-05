# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 14:52:01 2018

@author: swooty
"""

import numpy as np
import mpmath as mpm
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

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb = 6.6422e-26 #mass of potassium 40 in kg
k_B = 1.38064852e-23 # Boltzmann constant Joules per degree Kelvin

omega_x = 50.0*2.0*np.pi #Hz
omega_y = 50.0*2.0*np.pi #Hz
omega_z = 50.0*2.0*np.pi #Hz


N=1.0e6 #atom number 

omega_bar = (omega_x*omega_y*omega_z)**(1./3.)
ho_E = hbar*omega_bar

E_F = (6.0*N)**(1./3.) # in h.o. units 

def n_TF(E_ho,T_ho,zeta):
    n_E = zeta/(np.exp(E_ho/T_ho)+zeta)
    g_E = E_ho**2.0
    return n_E*g_E/2.0


zetaList = np.linspace(0.0,3000000000.0,num=5000)
polyLogList = np.array([float(mpm.polylog(3,-zeta)) for zeta in zetaList])

#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.5,hspace=1.0)
pan = fig.add_subplot(gs[0])
#pan.plot(zetaList,polyLogList,'b-')

#T=87.2e-9
#ThermE = T*k_B/ho_E #in h.o. units
#zeta = np.interp(N,-(ThermE**3.0)*polyLogList,zetaList)
#pan.plot(zetaList,np.abs(N+(ThermE**3.0)*polyLogList))
#print zeta

Tlist=np.array([436.0e-9,87.2e-9,21.8e-9]) #temperature in Kelvin 
zList = np.array([0.17016,30176.75,412346656.684])
label_list = np.array([r'$T = T_F$',r'$T = 0.2 T_F$',r'$T = 0.05 T_F$' ])

for i,T in enumerate(Tlist):
    ThermE = T*k_B/ho_E #in h.o. units
    zeta = zList[i]
    print zeta 
    
    bin_size = 1.0 # in h. o. units
    E_ho_list = bin_size*np.arange(600.0/bin_size)+bin_size/2
    
    pan.plot(E_ho_list,zeta/(np.exp(E_ho_list/ThermE)+zeta),label=label_list[i])
    #pan.plot(E_ho_list,n_TF(E_ho_list,ThermE,zeta)*bin_size/1e6,'r-',label='T=%.1f nK' %(T*1e9))
    print np.sum(n_TF(E_ho_list,ThermE,zeta)*bin_size)
    print zeta/(1.0+zeta)
    
nAtT0 = np.zeros(E_ho_list.size)
nAtT0[E_ho_list<E_F]= 1.0
pan.plot(E_ho_list,nAtT0,label=r'$T=0$')
pan.set_xlabel(r'Energy $\epsilon$ [$\hbar\bar{\omega}$]')
pan.set_ylabel(r'$n(\epsilon)$')
plt.legend()

plt.savefig('FermiDist.pdf',transparent=True)