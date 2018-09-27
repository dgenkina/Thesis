# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:04:11 2016

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

omega0=0.0
delta0=0.0
epsilon0=0.02
U0=5.0
k0=0.0
m0=1

omegaMax=0.5
rampOnt=0.038*Erecoil/hbar
tf=0.01*Erecoil/hbar
steps=600

def params(t):
    ramanRampOntr=0.018*Erecoil/hbar
    rampDetuning=0.02*Erecoil/hbar
    if t<ramanRampOntr:
        omega=omega0 + (omegaMax-omega0)*t/ramanRampOntr
        delta=delta0
        epsilon=epsilon0
        U=U0
    elif t<(ramanRampOntr+rampDetuning):
        omega=omegaMax
        delta=delta0*(1.0-(t-ramanRampOntr)/rampDetuning)
        epsilon=epsilon0
        U=U0
    else:
        omega=omegaMax
        delta=0.0
        epsilon=epsilon0
        U=U0
    
    return omega,delta,epsilon,U
    
def rampedOnPsi2(k0,omega0,delta0,epsilon0,U0,n,rampOnt,S,m0,steps=400):
    tlist=np.linspace(0.0,rampOnt,steps)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = LA.eig(RamanLatHam(k0, omega0, delta0, epsilon0, U0, n, S, m0))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    psi0=V1sorted[:,0]
    
    for tloc in tlist:
        omega,delta,epsilon,U=params(tloc)
        Energy,V=LA.eig(RamanLatHam(k0, omega, delta, epsilon, U, n, S, m0))
        V=V+0.0j
        Vinv=np.conj(V.transpose())
        psi0=np.dot(Vinv,psi0)
        teo=np.diag(np.exp(-1.0j*Energy*dt))
        psi0=np.dot(teo,psi0)
        psi0=np.dot(V,psi0)
        
  #  psi0=psi0.reshape(n,3).transpose()
   # pops=np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return psi0

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

def RamanLatHam(k, omega, delta, epsilon, U, n, S, m0):
    c=1064.0/790.0
    Nlat=2*n+1
    Ntot=Nlat*(2*S+1)
    Kinetic=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i%(2*S+1)-S)
        latI=np.float(np.divide(i,2*S+1)-n)
        Kinetic[i]=(k-2.0*(spinI-m0)*c-2.0*latI)**2.0
    H=np.diag(Kinetic)
    H+=delta*sLA.block_diag(*[Fz(S)]*Nlat)
    H+=((-1.0)**(S+1))*epsilon*sLA.block_diag(*[Fz(S)**2.0]*Nlat)
    H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
        
    for i in range(Ntot):
        for j in range(Ntot):
            if np.abs(i-j)==(2*S+1):
                H[i,j]=U/4.0         
    return H

def propagateRLHamiltonian2(t, k, omega, delta, epsilon, U, n,rampOnt,S,m0,**kwargs):  
    t=np.array(t)    
    psi0=rampedOnPsi2(k,omega,delta,epsilon,U,n,rampOnt,S,m0,**kwargs)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    omegaf,deltaf,epsilonf,Uf=params(rampOnt)
    H = RamanLatHam(k, omegaf, deltaf, epsilonf, Uf, n, S, m0)
    Energy, V = LA.eig(H)
    
    mat=np.array([np.identity(Energy.size)]*t.size)

    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    V=np.array([V]*t.size)
    Vinv=np.array([Vinv]*t.size)
    psi0=np.array([psi0]*t.size)


    U = np.einsum('ij,ijk->ijk',np.exp(-1j*np.outer(t, Energy)),mat)

    a = np.einsum('ijk,ik->ij',Vinv,psi0)
    b = np.einsum('ijk,ik->ij',U,a)                                                           
    psi = np.einsum('ijk,ik->ij',V,b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    latPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1),axis=1).flatten() #[:,0]#
    #populations in each spin state, summed over lattice sites, in time step blocks 
    return spinPops
    
def plotPulsedPops(tf, k, omega, delta, epsilon, U, n,rampOnt,S=2,m0=0,**kwargs):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    spinpops=propagateRLHamiltonian2(tlist, k, omega, delta, epsilon, U, n,rampOnt,S,m0,**kwargs)
    spinpops=spinpops.reshape(tlist.size,2*S+1).transpose()

    omegaf,deltaf,epsilonf,Uf=params(rampOnt)
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omegaf,deltaf,Uf))
    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[0],'m-',label=r'$m_F$=-2')
    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[1],'r-',label=r'$m_F$=-1')
    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[2],'g-',label=r'$m_F$=0')
    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[3],'b-',label=r'$m_F$=+1')
    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[4],'c-',label=r'$m_F$=+2')
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    return
    
plotPulsedPops(tf, k0, omega0, delta0, epsilon0, U0, 7,rampOnt,m0=m0,steps=800)