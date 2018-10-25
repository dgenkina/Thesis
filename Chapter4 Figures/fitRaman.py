# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 20:42:26 2016

@author: dng5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 
from numpy import linalg as LA
from matplotlib import cm
import readIgor


hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in nm
lambdaL = 1064e-9 # Lattice wavelength in nm
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaR**2.0) #recoil energy

S=1
m0=0
k=0.0
epsilon=0.0375


def Fz(S):
    a=np.arange(np.float(-S),np.float(S+1))
    F=np.diag(a)
    return F
    
def Fx(S):
    F=np.zeros((2*S+1,2*S+1))
    for i in range(int(2*S+1)):
        for j in range(int(2*S+1)):
            if np.abs(i-j)==1:
                F[i,j]=(1.0/2.0)*np.sqrt(S*(S+1)-(i-S)*(j-S))
    return F
    
def FxFlat(S):
    F=np.zeros((2*S+1,2*S+1))
    for i in range(2*S+1):
        for j in range(2*S+1):
            if np.abs(i-j)==1:
                F[i,j]=1.0/np.sqrt(2)
    return F 
    
def RamanF2(k, omega, delta, epsilon):
    Kinetic=(k+2.0*(np.arange(2*S+1)-S))**2.0
    H=np.diag(Kinetic)
    H+=delta*Fz(S)
    H+=((-1.0)**(S+1))*epsilon*(Fz(S)**2.0)
    H+=(np.sqrt(2.0)/2.0)*omega*Fx(S)   

    return H
    
def plotRamanBandStruct(omega, delta, epsilon,kList=np.linspace(-1.0,1.0,600),save=False,plot=True):
    i=0  
        

    s=2*S+1
    E=np.zeros((kList.size,s))

    pops=np.zeros((kList.size,s))
    m=np.zeros((kList.size,s))
    mFvect=np.array([np.float(j-S) for j in range(s)])

    for k in kList:
        H=RamanF2(k, omega, delta, epsilon)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()

        popsA=ev*np.conjugate(ev)
     #   pops[i]=popsA
        m[i]=np.dot(popsA,mFvect)

        i=i+1
    
    if plot:
        figure=plt.figure()
        panel=figure.add_subplot(1,1,1)
      #  panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        for i in range(s):   
            d=panel.scatter(kList,E[:,i],c=m[:,i],cmap=cm.jet_r,vmin=-S,vmax=S, marker='.')
        panel.set_xlabel(r'Momentum [$k_R$]',size=15)
        panel.set_ylabel(r'Energy [$E_R$]',size=15)
        cbar = figure.colorbar(d,ticks=np.arange(-S,S+1))
        cbar.ax.tick_params(labelsize=15) 
        cbar.set_label(r'$\langle m \rangle$',size=15)
      
        plt.tight_layout()
    return E,m
    
def propagateHamiltonian(t, psi0, k, omega, delta, epsilon):
    H=RamanF2(k,omega,delta,epsilon)
    #print H
    Energy, V = LA.eig(H)

    V=V+1j*0.0
 #   print V[:,0]

    Vinv = np.conjugate(np.transpose(V))
   # print Vinv
    #print np.dot(Vinv[1,:],(np.dot(H,V[:,1])))
    #print np.transpose(V)
    U=np.diag(np.exp(1j*np.array(Energy)*t))

    a=np.dot(Vinv,psi0)
   # print a
    b=np.dot(U,a)
    #print b
    psi = np.dot(V,b)
    #print psi
    pops=np.absolute(psi)**2.0

    return pops
    
def propagateHamiltonianFit(t,  omega, delta):
    t=np.array(t)
    psi0=np.zeros((t.size,2*S+1),dtype=complex)
    psi0[:,m0+S]=1.0+0.0j   


    H=RamanF2(k,omega,delta,epsilon)
    #print H
    Energy, V = LA.eig(H)


    mat=np.array([np.identity(Energy.size)]*t.size)
  
    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    

    V=np.array([V]*t.size)
    
    Vinv=np.array([Vinv]*t.size)


    U = np.einsum('ij,ijk->ijk',np.exp(-1j*np.outer(t, Energy)),mat)

    a = np.einsum('ijk,ik->ij',Vinv,psi0)
    b = np.einsum('ijk,ik->ij',U,a)                                                           
    psi = np.einsum('ijk,ik->ij',V,b)
    
    pops=np.absolute(psi)**2.0                     


    return pops.flatten()

if __name__=='__main__':
    roiList=np.array([[420,465,495,555],
            [530, 575, 585,645], 
             [615,660, 685,745]])
    
    filestart=39
    filestop=67
    fileroot = 'C:/Users/swooty/Documents/Thesis Data/Raman pulsing F=1/PIXIS_17Feb2016' 
    counts,fractions,waveDict,probe = readIgor.batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc="top")          
                
    tList=waveDict['holdDelay']
    tRecoils = tList*Erecoil/hbar
    fractions=np.array(fractions)

    popt,pcov = optimize.curve_fit(propagateHamiltonianFit,tRecoils,fractions.flatten(), p0=(4.0,0.0))

    tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),700)
    pops_fitted=propagateHamiltonianFit(tForFit,*popt).reshape(tForFit.size,2*S+1)

    
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title('Omega = ' + str(np.round(popt[0],3)) + ' Er/hbar, delta = '+str(np.round(popt[1],3))+' Er/hbar')
    panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,0],'bo', label='mF=+1')
    panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,1],'go', label='mF=0')
    panel.plot(tRecoils*hbar*1.0e6/Erecoil,fractions[:,2],'ro', label='mF=-1')
    panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[:,0],'b-')
    panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[:,1],'g-')
    panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[:,2],'r-')
    panel.set_xlabel('Raman hold time [us]')
    plt.legend()
    
    np.savez('RamanPulsingF1',time=tRecoils*hbar/Erecoil,tForFit=tForFit*hbar/Erecoil,
             fractions=fractions,pops_fitted=pops_fitted,popt=popt,pcov=np.sqrt(np.diag(pcov)))
    
