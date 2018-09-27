# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 12:23:41 2016

@author: dng5
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 
from numpy import linalg as LA

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in nm
lambdaL = 1064.0e-9 # Lattice wavelenght in nm
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
m0=-1

def RfHamiltonian(k, omega, delta, epsilon):
    H = np.array([[(k)**2.0 - delta, omega/2.0, 0.0],[omega/2.0, k**2.0-epsilon,omega/2.0],[0.0,omega/2.0,(k)**2.0 + delta]])
    return H
    
from scipy.linalg import block_diag

def propagateRfHamiltonian(t, omega, delta,epsilon):  
    t=np.array(t)
    k = 0.0
  #  epsilon=0.0
    psi0 = np.array([0+1j*0.0, 0.0+1j*0.0, 0.0+1j*0.0])
    psi0[m0+1]=1.0+0.0j
    H = RfHamiltonian(k, omega, delta ,epsilon)
    Energy, V = LA.eig(H)

    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))

    # np.outer(t, Energy).flatten() creates a matrix for all t
    U = np.diag(np.exp(-1j*np.outer(t, Energy).flatten()))  

    a = np.dot(Vinv, psi0)
    # This repeats a so that the shape is consitent with U
    aa = np.outer(np.ones(t.size),a).flatten()
                      
    # Have to add the transpose to make shapes match 
    b = np.dot(U, aa)                                     
    # Same block diagonal trick for eigenvector matrix
    VV = block_diag(*([V]*t.size))                          
    psi = np.dot(VV, b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    return pops,psi

  
def plotPulsedRf(tf, step, psi0 = [0.0,1.0,0.0], k=0.0, omega=4.0, delta=0.0, epsilon=0.0):
    tlist = np.arange(0,tf,step)
    pop0 = np.zeros(tlist.size)
    pop1 = np.zeros(tlist.size)
    pop2 = np.zeros(tlist.size)

    i=0    
    for t in tlist:
        pops,psi=propagateRfHamiltonian(t,omega,delta,epsilon)
     
        pop0[i]=pops[0]
        pop1[i]=pops[1]
        pop2[i]=pops[2]
        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.set_title('Omega = ' + str(omega) + ' E_L/hbar, delta = '+str(delta)+' E_L/hbar, epsilon = ' + str(epsilon)+ ' E_L')
    pan.plot(tlist*hbar/Erecoil,pop0,'r-', label='mF=-1')
    pan.plot(tlist*hbar/Erecoil,pop1,'g-', label='mF=0')
    pan.plot(tlist*hbar/Erecoil,pop2,'b-', label='mF=+1')
    pan.set_xlabel('Time [s]')
    legend()
    return 
    
#roiList=[[500, 550, 510,560], [500, 550, 460, 510] ,[500, 550,400,450]]
##
#filestart=73
#filestop=92
#fileroot = 'Y:/Data/2017/June/28/PIXIS_28Jun2017' 
##counts, fractions, waveDict, probeAvg =batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc='top')                      
#tList=waveDict['PulseTime']-0.000005
#tRecoils = tList*Erecoil/hbar
#fractions=np.array(fractions)
#
#a=np.array(tRecoils)
#
#popt,pcov = optimize.curve_fit(propagateRfHamiltonian,a,fractions.flatten(), p0=(4.2,0.0))
#print popt,pcov
#tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),200)
#pops_fitted=propagateRfHamiltonian(tForFit,*popt)
#sort=np.argsort(tRecoils)
#tSorted=tRecoils[sort]
#pop0 = np.array([pops_fitted[i*3] for i in range(tForFit.size)])
#pop1 = np.array([pops_fitted[i*3+1] for i in range(tForFit.size)]) 
#pop2 = np.array([pops_fitted[i*3+2] for i in range(tForFit.size)]) 
#
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.set_title(r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/\hbar$, $\delta$ = '+str(np.round(popt[1],3))+r' $E_L/\hbar$')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
#panel.plot(tList*1e6,fractions[:,0],'ro', label='mF=-1') #tRecoils*hbar*1e6/Erecoil
#panel.plot(tList*1e6,fractions[:,1],'go', label='mF=0')
#panel.plot(tList*1e6,fractions[:,2],'bo', label='mF=+1')
#panel.plot(tForFit*hbar*1e6/Erecoil,pop0,'r-')
#panel.plot(tForFit*hbar*1e6/Erecoil,pop1,'g-')
#panel.plot(tForFit*hbar*1e6/Erecoil,pop2,'b-')
#panel.set_xlabel(r'pulse time [$\mu s$]')
#legend()