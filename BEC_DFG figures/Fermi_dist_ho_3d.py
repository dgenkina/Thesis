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

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb = 6.6422e-26 #mass of potassium 40 in kg
k_B = 1.38064852e-23 # Boltzmann constant Joules per degree Kelvin

omega_x = 50.0*2.0*np.pi #Hz
omega_y = 50.0*2.0*np.pi #Hz
omega_z = 50.0*2.0*np.pi #Hz


N=1.0e6 #atom number 

omega_bar = (omega_x*omega_y*omega_z)**(1./3.)
ho_E = hbar*omega_bar



def n_TF(E_ho,T_ho,zeta):
    n_E = zeta/(np.exp(E_ho/T_ho)+zeta)
    g_E = E_ho**2.0
    return n_E*g_E/2.0


zetaList = np.linspace(0.999,1.0,num=5000)
polyLogList = np.array([float(mpm.polylog(3,-zeta)) for zeta in zetaList])

fig=plt.figure()
pan=fig.add_subplot(111)
pan.plot(zetaList,polyLogList,'b-')

T=206.0e-9
ThermE = T*k_B/ho_E #in h.o. units
zeta = np.interp(N,-(ThermE**3.0)*polyLogList,zetaList)
pan.plot(zetaList,np.abs(N+(ThermE**3.0)*polyLogList))
print zeta

#Tlist=np.array([436.0e-9,250.0e-9,40.0e-9]) #temperature in Kelvin 
#
#for T in Tlist:
#    ThermE = T*k_B/ho_E #in h.o. units
#    zeta = np.interp(N,-(ThermE**3.0)*polyLogList,zetaList)
#    print zeta 
#    
#    bin_size = 1.0 # in h. o. units
#    E_ho_list = bin_size*np.arange(3000.0/bin_size)+bin_size/2
#    
#    pan.plot(E_ho_list,zeta/(np.exp(E_ho_list/ThermE)+zeta),label='T=%.1f nK' %(T*1e9))
#    #pan.plot(E_ho_list,n_TF(E_ho_list,ThermE,zeta)*bin_size/1e6,'r-',label='T=%.1f nK' %(T*1e9))
#    print np.sum(n_TF(E_ho_list,ThermE,zeta)*bin_size)
#    print zeta/(1.0+zeta)
#plt.legend()