# -*- coding: utf-8 -*-
"""
Created on Wed May 24 12:38:11 2017

@author: dng5
"""
import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
c=4.0/3.0#0.0#1064.0/790.0#
Flat=True
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
    Nspin=int(2*S+1)
    F=np.zeros((Nspin,Nspin))
    for i in range(Nspin):
        for j in range(Nspin):
            if np.abs(i-j)==1:
                F[i,j]=1.0/np.sqrt(2)
    return F

def RamanLatHam(k, omega, delta, epsilon, U, n, S, m0):
    
    Nlat=2*n+1
    Nspin=int(2*S+1)
    Ntot=int(Nlat*(Nspin))
    Kinetic=np.zeros(Ntot)
    kStates=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i%(Nspin)-S)
        latI=np.float(np.divide(i,Nspin)-n)
        kStates[i]=k-2.0*(spinI-m0)*c-2.0*latI
        Kinetic[i]=(k-2.0*(spinI-m0)*c-2.0*latI)**2.0
    H=np.diag(Kinetic)
    H+=delta*sLA.block_diag(*[Fz(S)]*Nlat)
    try:
        H+=((-1.0)**(S+1))*epsilon*sLA.block_diag(*[Fz(S)**2.0]*Nlat)
    except ValueError:
        H+=epsilon*sLA.block_diag(*[Fz(S)**2.0]*Nlat)
    if Flat:
        H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[FxFlat(S)]*Nlat)    
    else:
         H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
    for i in range(Ntot):
        for j in range(Ntot):
            if np.abs(i-j)==(Nspin):
                H[i,j]=U/4.0         
    return H, kStates

n=11
S=0
#U=6.0
kList=np.linspace(-1.0,1.0,300)    
s=int(2*S+1)
N=2*n+1
E=np.zeros((kList.size,s*N))

mFvect=np.array([np.float(j-S) for j in range(s)])

Ulist = np.linspace(1.0,6.0,2)
J0list=np.zeros(Ulist.size)
J1list=np.zeros(Ulist.size)
tbGoodList=np.zeros(Ulist.size)
for ind,U in enumerate(Ulist):        
    i=0  
    for k in kList:
        H,kstates=RamanLatHam(k, 0.0, 0.0, 0.0, U, n,0,0)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        i=i+1
        
    #kgrid=np.fft.fftfreq(kList.size,d=(1.0)/(2.0*np.pi))    
    Jlist=np.fft.fft(E[:,0])/kList.size
    J0list[ind]=Jlist[0]
    J1list[ind]=Jlist[1]
    tbGoodList[ind]=(2.0*Jlist[1]**2.0)/(np.sum(np.dot(Jlist[1:],Jlist[1:])))

recover=np.zeros(kList.size)
inds=[0,1,-1]
for i in range(kList.size):
    recover[i]=np.sum(Jlist[inds].real*np.exp(np.arange(kList.size)[inds]*2.0*np.pi*i*1.0j/kList.size))
fig=plt.figure()
pan=fig.add_subplot(1,1,1)
pan.plot(kList,E[:,0],'b-')
pan.plot(kList,recover,'k-')
pan.plot(kList,Jlist[0]-2.0*Jlist[1]*np.cos(np.pi*kList),'m-')
fig=plt.figure()
pan=fig.add_subplot(1,1,1)
pan.plot(Ulist,J1list,'b-')
pan.set_xlabel(r'Lattice depth [$E_L$]',size='large')
pan.set_ylabel(r'Nearest neighbor tunneling [$E_L$]',size='large')

fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.plot(Ulist,tbGoodList,'b-')
pan2.set_xlabel(r'Lattice depth [$E_L$]',size='large')
pan2.set_ylabel(r'Fraction of band in 1st fourier component',size='large')

