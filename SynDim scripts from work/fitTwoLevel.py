# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 12:10:53 2017

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

def TwoLevelHamiltonian(k, omega, delta):
    H = np.array([[(k)**2.0 - delta/2, omega/2.0],[omega/2.0,(k)**2.0 + delta/2.0]])
    return H
    
from scipy.linalg import block_diag

def propagateTwoLevelHamiltonian(t, omega, delta):  
    k = 0.0
    psi0 = np.array([1.00+1j*0.0, 0.0+1j*0.0])
    H = TwoLevelHamiltonian(k, omega, delta)
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
    return pops
  
def plotPulsedTwoLevel(tf, step, psi0 = [1.0,0.0], k=0.0, omega=4.0, delta=0.0):
    tlist = np.arange(0,tf,step)
    pop0 = np.zeros(tlist.size)
    pop1 = np.zeros(tlist.size)

    i=0    
    for t in tlist:
        p0,p1=propagateTwoLevelHamiltonian(t,omega,delta)#,epsilon)
     
        pop0[i]=p0
        pop1[i]=p1

        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.set_title('Omega = ' + str(omega) + ' E_L/hbar, delta = '+str(delta)+' E_L/hbar')
    pan.plot(tlist*hbar/Erecoil,pop0,'b-', label='p0')
    pan.plot(tlist*hbar/Erecoil,pop1,'g-', label='p1')

    pan.set_xlabel('Time [s]')
    legend()
    return 

    
roiList=[[540, 580, 540, 580], [540, 580, 440, 480]]

filestart=34
filestop=63
fileroot = 'Y:/Data/2017/March/06/PIXIS_06Mar2017'  
counts, fractions, waveDict, probeAvg = batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc='top')                      
tList=waveDict['MicroFreqTrans1']
tRecoils = tList*Erecoil/hbar
fractions=np.array(fractions)

a=np.array(tRecoils)

#popt,pcov = optimize.curve_fit(propagateTwoLevelHamiltonian,a,fractions.flatten(), p0=(0.05,0.17))
#print popt,pcov
#tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),200)
#pops_fitted=propagateTwoLevelHamiltonian(tForFit,*popt)
#sort=np.argsort(tRecoils)
#tSorted=tRecoils[sort]
pop0 = np.array([pops_fitted[i*2] for i in range(tForFit.size)])
pop1 = np.array([pops_fitted[i*2+1] for i in range(tForFit.size)]) 


figure=plt.figure()
panel=figure.add_subplot(1,1,1)
#panel.set_title(r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/\hbar$, $\delta$ = '+str(np.round(popt[1],3))+r' $E_L/\hbar$')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
panel.plot(tList*1e6,fractions[:,0],'bo', label='|1,-1>') #tRecoils*hbar*1e6/Erecoil
panel.plot(tList*1e6,fractions[:,1],'go', label='|2,-1>')
#panel.plot(tForFit*hbar*1e6/Erecoil,pop0,'b-')
#panel.plot(tForFit*hbar*1e6/Erecoil,pop1,'g-')
panel.set_xlabel(r'Uwave frequency [MHz]')
legend()