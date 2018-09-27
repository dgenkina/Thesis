# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 18:06:48 2016

@author: dng5
"""
import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import stats
import matplotlib.pyplot as plt

S=2
m0=0


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
    Kinetic=(k-2.0*(np.arange(2*S+1)-S))**2.0
    H=np.diag(Kinetic)
    H+=delta*Fz(S)
    H+=((-1.0)**(S+1))*epsilon*(Fz(S)**2.0)
    H+=(np.sqrt(2.0)/2.0)*omega*Fx(S)   

    return H
    
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
    
kList = np.linspace(-2.0,2.0,300)
Energies=np.zeros((kList.size,2*S+1))
delta = 0.0
epsilon = 0.0
omega=1.0
i=0
for k in kList:
    H=RamanF2(k, omega, delta, epsilon)
    Energies[i], eigenstates = LA.eig(H)
    i=i+1

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
panel.plot(kList, Energies, 'b-')
  
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

def plotPulsedRaman(tf, step, k=0.0, omega=4.0, delta=0.0, epsilon=0.0):
    psi0=np.zeros(2*S+1,dtype=complex)
    psi0[m0+S]=1.0+0.0j    
    tlist = np.arange(0,tf,step)
    pops = np.zeros((tlist.size,2*S+1))

    i=0    
    for t in tlist:
        pops[i]=propagateHamiltonian(t, psi0,k,omega,delta,epsilon)
        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(tlist*hbar/Erecoil,pops[:,i],cDict[i-S], label='mF='+str(i-S))
    pan.set_title(r'$\Omega$=%.3f, $\delta$=%.3f,$\epsilon$=%.2f,k=%.2f'%(omega,delta,epsilon,k))
    legend()
    return 
    
def plotRamanOfDelta(deltaMin,deltaMax,step,t,k,omega,epsilon):
    psi0=np.zeros(2*S+1,dtype=complex)
    psi0[m0+S]=1.0+0.0j    
    dlist = np.arange(deltaMin,deltaMax,step)
    pops = np.zeros((dlist.size,2*S+1))
    trecoil=t*Erecoil/hbar

    i=0    
    for d in dlist:
        pops[i]=propagateHamiltonian(trecoil, psi0,k,omega,d,epsilon)
        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(dlist,pops[:,i],cDict[i-S], label='mF='+str(i-S))
    pan.set_title(r't=%.6f,$\Omega$=%.3f, $\epsilon$=%.2f,k=%.2f'%(t,omega,epsilon,k))
    pan.set_xlabel(r'Detuning [$E_R$]')
    legend()
    return 