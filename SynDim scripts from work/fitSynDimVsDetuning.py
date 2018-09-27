# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 18:44:46 2016

@author: dng5
"""

import numpy as np
from scipy import linalg as sLA
from scipy import optimize
from scipy import stats
import matplotlib.pyplot as plt

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
c=1064.0/790.0#
Flat=False
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
    if Flat:
        H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[FxFlat(S)]*Nlat)    
    else:
         H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
    for i in range(Ntot):
        for j in range(Ntot):
            if np.abs(i-j)==(2*S+1):
                H[i,j]=U/4.0         
    return H

def adiabaticPops(k, omega, delta, epsilon, U, n, S, m0):
    H=RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energies, eigenstates = sLA.eig(H)
    sort = np.argsort(Energies)
    Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
    ev1=eVsorted[:,0].reshape(2*n+1,2*S+1)
    m1pops=np.sum(ev1*np.conjugate(ev1),axis=0)#np.array([np.dot(ev1[:,0],ev1[:,0]),np.dot(ev1[:,1],ev1[:,1]),np.dot(ev1[:,2],ev1[:,2])])
    return m1pops
    
def plotAdiabaticPopsOfDelta(deltaMin,deltaMax,deltaStep,k,omega,epsilon,U,n,S,m0,rampOnt=0.02*Erecoil/hbar,holdt=0.0*Erecoil/hbar):
    dlist=np.arange(deltaMin,deltaMax,deltaStep)
    fracsA=np.zeros((dlist.size,2*S+1))
    fracsR=np.zeros((dlist.size,2*S+1))

    for ind,d in enumerate(dlist):
        fracsA[ind]=adiabaticPops(k, omega, d, epsilon, U, n, S, m0)
        fracsR[ind]=rampedOnPops(k,omega,d,epsilon,U,n,S,m0,rampOnt,holdt)

#    grad=np.gradient(fP)[dlist.size/2]/(dlist[1]-dlist[0])
#    print grad
#    print dlist[dlist.size/2]
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for mF in np.arange(-S,S+1):
        pan.plot(dlist,fracsA[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Adiabatic populations')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for mF in np.arange(-S,S+1):
        pan.plot(dlist,fracsR[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Pulsed populations')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    return dlist,fracsA,fracsR
    
def plotAdiabaticPopsOfOmega(omegaMin,omegaMax,omegaStep,k,delta,epsilon,U,n,S,m0,rampOnt=0.02*Erecoil/hbar,holdt=0.0*Erecoil/hbar):
    olist=np.arange(omegaMin,omegaMax,omegaStep)
    fracsA=np.zeros((olist.size,2*S+1))
    fracsR=np.zeros((olist.size,2*S+1))

    for ind,o in enumerate(olist):
        fracsA[ind]=adiabaticPops(k, o, delta, epsilon, U, n, S, m0)
        fracsR[ind]=rampedOnPops(k,o,delta,epsilon,U,n,S,m0,rampOnt,holdt)

#    grad=np.gradient(fP)[dlist.size/2]/(dlist[1]-dlist[0])
#    print grad
#    print dlist[dlist.size/2]
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for mF in np.arange(-S,S+1):
        pan.plot(olist,fracsA[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
    pan.set_xlabel(r'$\Omega$ [$E_r/h$]')
    pan.set_ylabel('Adiabatic populations')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for mF in np.arange(-S,S+1):
        pan.plot(olist,fracsR[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
    pan.set_xlabel(r'$\Omega$ [$E_r/h$]')
    pan.set_ylabel('Pulsed populations')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    return olist,fracsA,fracsR
    
def rampedOnPops(k,omega,delta,epsilon,U,n,S,m0,rampOnt,holdt):
    tlist=np.linspace(0.0,rampOnt,80)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = sLA.eig(RamanLatHam(k, 0.0, 0.0, 0.0, U, n, S, m0))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    psi0=V1sorted[:,0]

    
    
    for t in tlist:
        omegaLoc=omega*t/rampOnt
        Energy,V=sLA.eig(RamanLatHam(k,omegaLoc,delta,epsilon, U, n, S, m0))
        V=V+0.0j
        Vinv=np.conj(V.transpose())
        psi0=np.dot(Vinv,psi0)
        teo=np.diag(np.exp(-1.0j*Energy*dt))
        psi0=np.dot(teo,psi0)
        psi0=np.dot(V,psi0)
    
    Energy,V=sLA.eig(RamanLatHam(k,omega,delta,epsilon, U, n, S, m0))
    V=V+0.0j
    Vinv=np.conj(V.transpose())
    psi0=np.dot(Vinv,psi0)
    teo=np.diag(np.exp(-1.0j*Energy*holdt))
    psi0=np.dot(teo,psi0)
    psi0=np.dot(V,psi0)    
    
        
    psi0=psi0.reshape(2*n+1,2*S+1)
    pops=np.sum(psi0*np.conjugate(psi0),axis=0)#np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return pops
    
def popsForFit(freqList,resFreq,omega,U):
    epsilon=0.048
    n=11
    k=0.0
 #   U=6.5
    wP=1.0
    wM=1.0
    weights=np.array([wP,1.0,wM])
    deltaList=(resFreq-freqList)
    rampOnt=0.08*Erecoil/hbar
    popsW=np.zeros([deltaList.size,3])
    for ind,delta in enumerate(deltaList):
        pop=rampedOnPops(k,omega,delta,epsilon,U,n,rampOnt)
        totW=np.dot(pop,weights)
        popsW[ind]=pop*weights/totW
        
    return popsW.flatten()
    
def popsForFitPenalized(freqList,resFreq,omega,wP,wM,U):
    popsW=popsForFit(freqList,resFreq,omega,wP,wM,U)
    if ((wP<1.3) & (wP>0.7)):
        pen1=0
    else:
        pen1=1.0
    if ((wM<1.3) & (wM>0.7)):
        pen2=0
    else:
        pen2=1.0
    if ((U<8.0)&(U>4.5)):
        pen3=0.0
    else:
        pen3=1.0
    return popsW-pen1-pen2-pen3

    
    
filename='18Jul2016_Detuning2.npz'
datafile=np.load(filename)
freqList=datafile['tlist'][2:14]
freqRecoils = freqList*1.0e6*(hbar*2.0*pi)/Erecoil
rfResGuess=0.81721*1.0e6*(hbar*2.0*pi)/Erecoil
fractionP=datafile['fractionP'][2:14]
fraction0=datafile['fraction0'][2:14]
fractionM=datafile['fractionM'][2:14]
fractions=np.array([fractionP,fraction0,fractionM]).transpose()

popt,pcov = optimize.curve_fit(popsForFit,freqRecoils,fractions.flatten(), p0=(rfResGuess,0.5,6.0))
print popt,pcov
freqForFit=np.linspace(np.min(freqRecoils),np.max(freqRecoils),60)
pops_fitted=popsForFit(freqForFit,*popt).reshape(freqForFit.size,3).transpose()
pop0 =pops_fitted[0]
pop1 = pops_fitted[1]
pop2 = pops_fitted[2]



figure=plt.figure()
panel=figure.add_subplot(1,1,1)
panel.set_title(r'$\Omega$ = ' + str(np.round(popt[1],2)) + ' $E_L/\hbar$, U ='+str(np.round(popt[2],3))+', resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(hbar*2.0*np.pi),5))+' MHz')
panel.plot(freqList,fractionP,'bo', label='mF=+1')
panel.plot(freqList,fraction0,'go', label='mF=0')
panel.plot(freqList,fractionM,'ro', label='mF=-1')
panel.plot(freqForFit*Erecoil*1e-6/(hbar*2.0*np.pi),pop0,'b-')
panel.plot(freqForFit*Erecoil*1e-6/(hbar*2.0*np.pi),pop1,'g-')
panel.plot(freqForFit*Erecoil*1e-6/(hbar*2.0*np.pi),pop2,'r-')
panel.set_xlabel('Raman difference frequency [MHz]')