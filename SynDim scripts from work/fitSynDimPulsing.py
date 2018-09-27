# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 17:42:37 2016

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy


def RamanLatHamiltonian(k, omega, delta, epsilon, U, n):
    n=np.int(n)
    if n%2==0:
        print "Number of lattice orders n must be odd!"
        
    c=1.064/0.79
    H=np.zeros((3*n,3*n))
    for i in np.arange(3*n):
        spinI=i%3-1
        latI=np.divide(i,3)-np.divide(n,2)
        for j in np.arange(3*n):
            spinJ=j%3-1
            if i==j:
                spinI=np.float(spinI)
                latI=np.float(latI)
                H[i,j]=(k-2.0*spinI*c-2.0*latI)**2.0+spinI*delta-(1-np.abs(spinI))*epsilon
            if np.abs(i-j)==3:
                H[i,j]=U/4.0
            if ((np.abs(i-j)==1) and (np.abs(spinI-spinJ)==1)):
                H[i,j]=omega/2.0
            
    return H
    
def propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n):  
    #k = 0.0
    #epsilon=0.0
    Energy1, V1 = LA.eig(RamanLatHamiltonian(k, 0.0, 0.0, epsilon, U, n))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    psi0=V1sorted[:,0]
 #   psi0[np.divide(3*n,2)]=1.0+0.0*1j
 #   psi0=np.zeros(3*n, dtype=complex)
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHamiltonian(k, omega, delta ,epsilon,U,n)
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
    VV = sLA.block_diag(*([V]*t.size))                          
    psi = np.dot(VV, b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    latPops=np.sum(pops.reshape(t.size,n,3)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,n,3),axis=1).flatten() 
    #populations in each spin state, summed over lattice sites, in time step blocks 
    return spinPops
    
roiList=np.array([[590, 620, 540, 585], 
         [520, 550, 535, 585],
         [450,480, 540, 585]])

filestart=79
filestop=108
fileroot = 'X:/2016/June/20/PIXIS_20Jun2016'  
#counts,fractions,waveDict = readIgor.batchCountMultipleROIs(fileroot,filestart,filestop,roiList)                      
#tList=waveDict['holdLat']
#tRecoils = np.array(tList*Erecoil/hbar)
#fractions=np.array(fractions)
#frac0=np.array(fractions[:,0])
#print tRecoils
#psi0=np.array([0+1j*0.0,1.0+1j*0.0,0.0+1j*0.0])
n=9
c=4.0/3.0
k=0.1
omega=0.65
delta=0.0
epsilon=0.048
U=6.0
#def constrainedPropagateHam( k, omega, delta, epsilon, n):
#    return lambda t, U: np.array(propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n))
#HofT=constrainedPropagateHam( k, omega, delta, epsilon, n)
##print H(1.0,4.0,0.0,0.0)
#a=np.array(tRecoils)
#b=np.array(frac0)
#popt,pcov = optimize.curve_fit(HofT,tRecoils,fractions.flatten(), p0=(12.0))
#print popt, pcov
tForFit=np.linspace(0.0,0.003*Erecoil/hbar,80)#np.linspace(np.min(tRecoils),np.max(tRecoils),40)
pops_fitted=propagateRLHamiltonian(tForFit,*(k, omega, delta, epsilon, U, n))
#sort=np.argsort(tRecoils)
#tSorted=tRecoils[sort]
pop1 = np.array([pops_fitted[i*3+2] for i in range(tForFit.size)])
pop0 = np.array([pops_fitted[i*3+1] for i in range(tForFit.size)]) 
pop2 = np.array([pops_fitted[i*3] for i in range(tForFit.size)]) 

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
panel.set_title(r'k=%.1f $k_L$,$\Omega$=%.1f $E_L$, $\delta$=%.1f $E_L$, $\epsilon$ = %.3f $E_L$ U =%.1f $E_L$'%(k, omega, delta, epsilon, U))
#panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,0],'bo', label='+2k_L')
#panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,1],'go', label='0')
#panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,2],'ro', label='-2K_L')
panel.plot(tForFit*hbar*1.0e3/Erecoil,pop0,'g-')
panel.plot(tForFit*hbar*1.0e3/Erecoil,pop1,'b-')
panel.plot(tForFit*hbar*1.0e3/Erecoil,pop2,'r-')
panel.set_xlabel('Lattice pulse time [ms]')
panel.set_ylabel('Fractional Populations')
legend()