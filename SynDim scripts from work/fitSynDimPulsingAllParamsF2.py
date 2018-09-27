# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:18:23 2016

@author: dng5
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:08:06 2016

@author: dng5
"""



import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy


c=1064.0/790.0#0.0#
n=7
k=0.0
epsilon=0.0207
U=4.4
S=1
m0=0
#delta=0.0
#omega=0.5
tau=0.0003*Erecoil/hbar
rampOnt=0.0001*Erecoil/hbar


#dataFile=np.load('07Mar2017_files_694-723.npz')
omegaGuess=0.5
deltaGuess=-0.03

cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

weight=False
Flat=False
w0=1.08
wP=1.12


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
    kStates=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i%(2*S+1)-S)
        latI=np.float(np.divide(i,2*S+1)-n)
        kStates[i]=k-2.0*(spinI-m0)*c-2.0*latI
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

def kStates(k, omega, delta, epsilon, U, n, S, m0):

    Nlat=2*n+1
    Ntot=Nlat*(2*S+1)

    kstates=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i%(2*S+1)-S)
        latI=np.float(np.divide(i,2*S+1)-n)
        kstates[i]=k-2.0*(spinI-m0)*c-2.0*latI
        
    return kstates

def expRamp(t,tMax,start,stop,tau):
    out=(stop-start*np.exp(tMax/tau)+(start-stop)*np.exp(t/tau))/(1-np.exp(tMax/tau))
    return out
    
def rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,rampOnt=0.0003*Erecoil/hbar,steps=50):
    tlist=np.linspace(0.0,rampOnt,steps)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = LA.eig(RamanLatHam(k, 0.0, 0.0, epsilon, U, n, S, m0))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    if c==0.0:
        psiLat=V1sorted[:,0].reshape(2*n+1,2*S+1).transpose()
        psi0=np.zeros(Energy1.size).reshape(2*n+1,2*S+1).transpose()
        psi0[m0+S]=psiLat[1]/np.sqrt(np.sum(psiLat[1]**2.0))
        psi0=psi0.transpose().flatten()
       
    else:
        psi0=V1sorted[:,0]
    
    for t in tlist:
        omegaLoc=expRamp(t,rampOnt,0.0,omega,tau)
        Energy,V=sLA.eigh(RamanLatHam(k, omegaLoc, delta, epsilon, U, n, S, m0))
        V=V+0.0j
        Vinv=np.conj(V.transpose())
        psi0=np.dot(Vinv,psi0)
        teo=np.diag(np.exp(-1.0j*Energy*dt))
        psi0=np.dot(teo,psi0)
        psi0=np.dot(V,psi0)
        
  #  psi0=psi0.reshape(n,3).transpose()
   # pops=np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return psi0

    
def propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n,S,m0):  
    t=np.array(t)    
    psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
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
    latPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1),axis=1).flatten() 
    #populations in each spin state, summed over lattice sites, in time step blocks 
    return spinPops
    
def propagateRLHamiltonian2(t, k, omega, delta, epsilon, U, n,S,m0,**kwargs):  
    t=np.array(t)    
    psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,**kwargs)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energy, V = LA.eig(H)
    
    mat=np.array([np.identity(Energy.size)]*t.size)


    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    V=np.array([V]*t.size)
    Vinv=np.array([Vinv]*t.size)
    psi0=np.array([psi0]*t.size)

    aa=np.exp(-1j*np.outer(t, Energy))

    U = np.einsum('ij,ijk->ijk',aa,mat)

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

def propagateRLHamiltonianAllPops(t, k, omega, delta, epsilon, U, n,S,m0,**kwargs):  
    t=np.array(t)    
    psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,**kwargs)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
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

    return pops
    
def plotPulsedPops(tf, k, omega, delta, epsilon, U, n,S=2,m0=0,**kwargs):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    spinpops=propagateRLHamiltonian2(tlist, k, omega, delta, epsilon, U, n,S,m0,**kwargs)
    spinpops=spinpops.reshape(tlist.size,2*S+1).transpose()

    
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omega,delta,U))
    for s in range(2*S+1):
        panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[s],label=r'$m_F$='+str(int(s-S)))
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[1],'r-',label=r'$m_F$=-1')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[2],'g-',label=r'$m_F$=0')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[3],'b-',label=r'$m_F$=+1')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[4],'c-',label=r'$m_F$=+2')
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    return
    
    
def plotPulsedPopsGen(tf, k, omega, delta, epsilon, U, n,orders=[(0,1),(-1,1),(1,1)],S=2,m0=0,**kwargs):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    popsAllraw=propagateRLHamiltonianAllPops(tlist, k, omega, delta, epsilon, U, n,S,m0,**kwargs)
    popsAll=popsAllraw.reshape(tlist.size,2*n+1,2*S+1)
    pops=np.zeros((tlist.size,len(orders)))
    for i in range(len(orders)):
        nn=n-orders[i][0]
        ss=S+orders[i][1]
        pops[:,i]=popsAll[:,nn,ss]

    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omega,delta,U))
    for i in range(len(orders)):
        panel.plot(tlist*hbar*1.0e6/Erecoil,pops[:,i],label=r'($k_L$,$m_F$)='+str(orders[i]))
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    
    #Plot normal spin populations 
    spinpops=np.sum(popsAllraw.reshape(tlist.size,2*n+1,2*S+1),axis=1).flatten()
    spinpops=spinpops.reshape(tlist.size,2*S+1).transpose()
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omega,delta,U))
    for s in range(2*S+1):
        panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[s],label=r'$m_F$='+str(int(s-S)))

    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    
    #Calcualte chiral current
    
    mFlist = np.arange(2*S+1) - S
    kstates = kStates(k, omega, delta, epsilon, U, n,S,m0)
    current = np.sum(np.dot((kstates*popsAllraw).reshape(tlist.size,2*n+1,2*S+1),mFlist),axis=1)
    
    figure = plt.figure()
    panel = figure.add_subplot(111)
    panel.plot(tlist*hbar*1.0e3/Erecoil,current)
    panel.set_xlabel('Lattice pulse time [ms]')
    panel.set_ylabel('Chiral current')
    
    
    return
    
def adiabaticPops(k,omega,delta,epsilon,U,n,S,m0):
    H=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
    Energies, eigenstates = LA.eig(H)
    sort = np.argsort(Energies)
    Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
    ev1=eVsorted[:,0].reshape(2*n+1,2*S+1)
    m1pops=np.einsum('ij,ij->j',ev1,np.conjugate(ev1))#np.array([np.dot(ev1[:,0],ev1[:,0]),np.dot(ev1[:,1],ev1[:,1]),np.dot(ev1[:,2],ev1[:,2])])
    return m1pops
    
def plotAdiabatiPopsOfDelta(deltaMin,deltaMax,deltaStep,k,omega,epsilon,U,n,S,m0,t,rampOnt=0.015*Erecoil/hbar):
    dlist=np.arange(deltaMin,deltaMax,deltaStep)
    fAd=np.zeros((dlist.size,2*S+1))
    fRamped=np.zeros((dlist.size,2*S+1))


    for ind,d in enumerate(dlist):
        pops=adiabaticPops(k,omega,d,epsilon,U,n,S,m0)
        fAd[ind]=pops
        pops2=propagateRLHamiltonian2(t, k, omega, d, epsilon, U, n,S,m0,rampOnt=rampOnt)
        fRamped[ind]=pops2


    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(dlist,fAd[:,i], cDict[i-S], label='mF='+str(i-S))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Fractional pop, adiabatic')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(dlist,fRamped[:,i],cDict[i-S], label='mF='+str(i-S))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Fractional pop, ramped')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,pulse time=%.3f ms,$\phi=$%.2f'%(omega,U,epsilon,k,n,t*hbar*1.0e3/Erecoil,c))
    legend()
    return

    
def plotAdiabatiPopsOfU(UMin,UMax,UStep,k,omega,epsilon,delta,n,S,m0,t,rampOnt=0.015*Erecoil/hbar):
    Ulist=np.arange(UMin,UMax,UStep)
    fAd=np.zeros((Ulist.size,2*S+1))
    fRamped=np.zeros((Ulist.size,2*S+1))


    for ind,U in enumerate(Ulist):
        pops=adiabaticPops(k,omega,delta,epsilon,U,n,S,m0)
        fAd[ind]=pops
        pops2=propagateRLHamiltonian2(t, k, omega, delta, epsilon, U, n,S,m0,rampOnt=rampOnt)
        fRamped[ind]=pops2


    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(Ulist,fAd[:,i], cDict[i-S], label='mF='+str(i-S))
    pan.set_xlabel(r'Lattice Depth [$E_r$]')
    pan.set_ylabel('Fractional pop, adiabatic')
    pan.set_title(r'$\Omega$=%.3f, $\delta$=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,delta,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(Ulist,fRamped[:,i],cDict[i-S], label='mF='+str(i-S))
    pan.set_xlabel(r'Lattice Depth [$E_r$]')
    pan.set_ylabel('Fractional pop, ramped')
    pan.set_title(r'$\Omega$=%.3f, $\delta$=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,pulse time=%.3f ms,$\phi=$%.2f'%(omega,delta,epsilon,k,n,t*hbar*1.0e3/Erecoil,c))
    legend()
    return
    
    
#signalGood=dataFile['signalGood']
#imbal=dataFile['imbalArray']
#cutoff=0.25
#fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
##fieldGoodArray=np.ones(dataFile['tlist'].size,dtype=bool)
#toffs=0.0
#tList=dataFile['tlist'][fieldGoodArray]+toffs
#sort=np.argsort(tList)
#tList=tList[sort]
#
#
#fractionP=dataFile['fractionP'][fieldGoodArray][sort]
#fraction0=dataFile['fraction0'][fieldGoodArray][sort]
#fractionM=dataFile['fractionM'][fieldGoodArray][sort]
#if weight:
#    tot=fractionM+w0*fraction0+wP*fractionP
#    fraction0=fraction0*w0/tot
#    fractionP=fractionP*wP/tot
#    fractionM=fractionM/tot
#
#maxInd=np.int(fractionP.size)
#fractions=np.append(np.append(fractionM[:maxInd],fraction0[:maxInd]),fractionP[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']   
#fractions=fractions.reshape(3,fraction0[:maxInd].size).transpose()
#
#if S==2:
#    fractionP2=dataFile['fractionP2'][fieldGoodArray][sort]
#    fractionM2=dataFile['fractionM2'][fieldGoodArray][sort]
#    fractions=np.append(np.append(np.append(np.append(fractionM2[:maxInd],fractionM[:maxInd]),fraction0[:maxInd]),fractionP[:maxInd]),fractionP2[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']   
#    fractions=fractions.reshape(5,fraction0[:maxInd].size).transpose()
#
#
##tList=waveDict['pulseDelay']#np.load('SynDimCalibData31Mar16.npz')['tList']
#tRecoils = np.array(tList*Erecoil/hbar)
#fractions=np.array(fractions)
#frac0=np.array(fractions[:,0])
##print tRecoils
##psi0=np.array([0+1j*0.0,1.0+1j*0.0,0.0+1j*0.0])
#
#
#def constrainedPropagateHam( k,epsilon,U, n,S,m0):
#    return lambda t,omega,delta: np.array(propagateRLHamiltonian2(t, k, omega, delta, epsilon, U,n,S,m0,rampOnt=rampOnt))
#HofT=constrainedPropagateHam( k,epsilon,U, n,S,m0)
###print H(1.0,4.0,0.0,0.0)
##a=np.array(tRecoils)
##b=np.array(frac0)
#
#popt,pcov = optimize.curve_fit(HofT,tRecoils[:maxInd],fractions.flatten(), p0=(omegaGuess,deltaGuess))
#print popt, np.sqrt(np.diag(pcov))
##popt=(0.6,0.08)
#tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),100)
#pops_fitted=HofT(tForFit,*popt)
##sort=np.argsort(tRecoils)
##tSorted=tRecoils[sort]
#
#pops_fitted=pops_fitted.reshape(tForFit.size,2*S+1).transpose()
#
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.set_title( r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/h$,$\delta$ = '+str(np.round(popt[1],3))+r' $E_L/h$')#, U= '+str(np.round(popt[2],3))+r'$E_L$')
#
#panel.plot(tList*1.0e6,fractionP,'bo', label=r'$m_F$=+1')
#panel.plot(tList*1.0e6,fraction0,'go', label=r'$m_F$=0')
#panel.plot(tList*1.0e6,fractionM,'ro', label=r'$m_F$=-1')
#
#
#panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[S],'g-')
#panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[S+1],'b-')
#panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[S-1],'r-')
#
#if S==2:
#    panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[S+2],'c-')
#    panel.plot(tForFit*hbar*1.0e6/Erecoil,pops_fitted[S-2],'m-')
#    panel.plot(tList*1.0e6,fractionM2,'mo', label=r'$m_F$=-2')
#    panel.plot(tList*1.0e6,fractionP2,'co', label=r'$m_F$=+2')
#panel.set_xlabel('Lattice pulse time [us]')
#legend()
##
##figure=plt.figure()
##pan=figure.add_subplot(1,1,1)
##t=np.linspace(0.0,rampOnt,50) 
##pan.plot(t ,expRamp(t,rampOnt,0.0,omega,tau),'b-')
#
##figure2=plt.figure()
##panel2=figure2.add_subplot(1,1,1)
##
##panel2.plot(imbal,fractionP,'bo', label=r'$m_F$=+1')
##panel2.plot(imbal,fraction0,'go', label=r'$m_F$=0')
##panel2.plot(imbal,fractionM,'ro', label=r'$m_F$=-1')
##
##panel2.set_xlabel('uwave imbalance')
##legend()
##print 'Fraction in +1 = '+ np.str(np.round(np.average(fractionP),4))+' +/-' + np.str(np.round(np.std(fractionP),4))
##print 'Fraction in 0 = '+ np.str(np.round(np.average(fraction0),4))+' +/-' + np.str(np.round(np.std(fraction0),4))
##print 'Fraction in -1 = '+ np.str(np.round(np.average(fractionM),4))+' +/-' + np.str(np.round(np.std(fractionM),4))