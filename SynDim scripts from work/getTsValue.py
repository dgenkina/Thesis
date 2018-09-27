# -*- coding: utf-8 -*-
"""
Created on Mon Jul 09 12:37:37 2018

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt

import tightBindingHam
import SynDimBandStructureGeneral


hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=7.0
cLoc=0.0#3.0/q#1064.0/790.0#
Flat=False

omega = 0.0
delta = 0.0
epsilon = 0.02
U = 4.4
n = 7
S = 1
s=2*S+1
kList,E,pops, mFull = SynDimBandStructureGeneral.plotSynDimBandStructGen(omega, delta, epsilon, U, n, S, 0,c=cLoc,save=False,magCell=False,plot=False)
Efull = E[:,0:5]
offset1 = np.average(E[:,0])

tx = 0.0775
Flat=True
epsilon=0.0
ts=0.1

tsList = np.linspace(0.45,0.55,10)
diffList = np.zeros(tsList.size)

for i,ts in enumerate(tsList):
    kList, Etb, mTB = tightBindingHam.plotSynDimBandStructGenTBClean(ts, delta, epsilon, tx, S,0,c=cLoc,kList=kList,save=False,plots=False)
    offset2=np.average(Etb[:,0])
    diffList[i]=np.sum(np.abs(Etb[:,0:2]-offset2+offset1-Efull[:,0:2]))
    
fig=plt.figure()
pan=fig.add_subplot(111)
pan.plot(tsList,diffList)

print np.min(diffList)
ts = tsList[np.argmin(diffList)]
print ts
kList, Etb, mTB = tightBindingHam.plotSynDimBandStructGenTBClean(ts, delta, epsilon, tx, S,0,c=cLoc,kList=kList,save=False,plots=False)
offset2=np.average(Etb[:,0])

figure=plt.figure()
panel=figure.add_subplot(111)

lw=5
d=panel.plot(kList/2.0,(Etb[:,0]-offset2+offset1)/tx,'b-',lw=lw,label = 'Tight binding') #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(kList/2.0,Efull[:,0]/tx,'r-',lw=lw, label = 'Full Hamiltonian') #c=mFull[:,i],vmin=-S,vmax=S,
d=panel.plot(kList/2.0,(Etb[:,1]-offset2+offset1)/tx,'b-',lw=lw) #c=mTB[:,i],vmin=-S,vmax=S, 
d=panel.plot(kList/2.0,Efull[:,1]/tx,'r-',lw=lw) #c=mFull[:,i],vmin=-S,vmax=S,
panel.set_xlabel(r'Crystal momentum $q_x$ [$2k_L$]')   
panel.set_ylabel(r'Energy [$t_x$]')
legend()
#plt.colorbar(d)
panel.set_title('tx = %.3f, ts = %.3f tx' %(tx,ts/tx))
    
panel.set_xlim(-0.5,0.5)