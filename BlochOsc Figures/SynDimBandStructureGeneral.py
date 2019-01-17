# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:34:49 2016

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=7.0
c=1064.0/790.0#0.0#3.0/q#
Flat=False


cdictC = {'red':   [(0.0,  0.0, 1.0), 
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)],

         'blue':  [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)]}
                   
cyans = mpl.colors.LinearSegmentedColormap('Cyans', cdictC)
plt.register_cmap(cmap=cyans)        

cdictM = {'red':   [(0.0,  0.0, 1.0), 
                   (1.0,  1.0, 0.0)],

         'green': [(0.0,  0.0, 1.0),
                   (1.0,  0.0, 1.0)],

         'blue':  [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)]}
                   
magentas = mpl.colors.LinearSegmentedColormap('Magentas', cdictM)
plt.register_cmap(cmap=magentas)



cmapDict={}
cmapDict[-2]='Magentas'
cmapDict[-1]='Reds'
cmapDict[0]='Greens'
cmapDict[1]='Blues'
cmapDict[2]='Cyans'


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

def RamanLatHam(k, omega, delta, epsilon, U, n, S, m0,c=c, Flat=Flat):
    
    Nlat=2*n+1
    Nspin=int(2*S+1)
    Ntot=int(Nlat*(Nspin))
    Kinetic=np.zeros(Ntot)
    kStates=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i%(Nspin)-S)
        latI=np.float(np.divide(i,Nspin)-n)
        kStates[i]=k+2.0*(spinI-m0)*c-2.0*latI
        Kinetic[i]=(k+2.0*(spinI-m0)*c-2.0*latI)**2.0
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
    
def plotSynDimBandStructGen(omega, delta, epsilon, U, n,S, m0,c=c,kList=np.linspace(-1.0,1.0,600),save=False,magCell=False,plot=True,Flat=Flat):
    i=0  
        
    s=int(2*S+1)
    N=2*n+1
    E=np.zeros((kList.size,s*N))
    m=np.zeros((kList.size,s*N))
    pops=np.zeros((kList.size,s*N,s))
    ktot=np.zeros((kList.size,s*N))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
    if magCell:
        magInUnitCell=np.zeros((kList.size,s-2))
    for k in kList:
        H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0,c=c,Flat=Flat)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(N*s,N,s)
        popsA=np.sum(evA*np.conjugate(evA),axis=1)
        pops[i]=popsA
        
        evB=ev.reshape(N*s,N*s)
        popsB=evB*np.conjugate(evB)
        ktot[i]=np.einsum('i,ji->j',kstates,popsB)
        
        m[i]=np.dot(popsA,mFvect)
        if magCell:
            for j in range(s-2):
                magInUnitCell[i,j]=(-popsA[0,j]+popsA[0,j+2])/(popsA[0,j]+popsA[0,j+1]+popsA[0,j+2])
        i=i+1
    
    if plot:
        figure=plt.figure()
        panel=figure.add_subplot(111)
        for i in range(s):   
            d=panel.scatter(kList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='.',lw=0)
        panel.set_xlim(-1.0,1.0)
#    fig2=plt.figure()
#    pan2=fig2.add_subplot(1,1,1)    
#    pan2.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
#
#    for i in range(s):   
#        pan2.plot(kList,ktot[:,i],label='band '+str(i))
#    pan2.set_xlabel(r'$q/k_L$')
#    pan2.set_ylabel(r'Total momentum $[k_L]$')
#    plt.legend()
    
    
#    kM=kList[np.where(pops[:,0,S-1]==np.max(pops[:,0,S-1]))]
#    k0=kList[np.where(pops[:,0,S]==np.max(pops[:,0,S]))]
#    kP=kList[np.where(pops[:,0,S+1]==np.max(pops[:,0,S+1]))]
#    print kM, k0, kP
#    
#    chern=2.0*2.0/(3.0*(kP-kM))
#    print 'chern number measured = ' +str(chern)
    
#    if plot: 
#        fig3=plt.figure()
#        #fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
#        pan3=fig3.add_subplot(1,1,1)    
#       # pan3.set_title('Lowest band to second band')    
#        for i in range(s):    
#            pan3.scatter(kList,[i-S for j in range(kList.size)],c=pops[:,0,i],cmap=cmapDict[S-i],vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
#       #     pan3.scatter(kList+2.0,[i for j in range(kList.size)],c=pops[:,1,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
#        
#        pan3.set_ylabel(r'$m$ state', size=20)
#        pan3.set_xlabel(r'Crystal momentum [$k_L$]',size=20)
#        plt.savefig('Z:/My Documents/papers etc/coloredBarcodesTheory2_chern2.jpg')

    #    
#    fig3=plt.figure()
#    fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
#    pan3=fig3.add_subplot(1,1,1)    
#    pan3.set_title('Lowest bands')  
#    
    popAvg=np.zeros(pops[:,0,:].shape)
    for ind,k in enumerate(kList):
        for i in range(s):
            popAvg[ind,i]=np.average(pops[ind,0:int(np.ceil(s/q))-1,i])
#        
#    for i in range(s):    
#        pan3.scatter(kList,[i for j in range(kList.size)],c=popAvg[:,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
#    
#    pan3.set_ylabel('Synthetic lattice site')
#    pan3.set_xlabel(r'Crystal momentum [$k_L$]')
    
#    maxSite=np.zeros(kList.size)
#    for i in np.arange(kList.size):
#        maxSite[i]=np.where(pops[i,0,:]==np.max(pops[i,0,:]))[0][0]
#        
#    fig4=plt.figure()
#    pan4=fig4.add_subplot(1,1,1)
#    pan4.plot(kList,maxSite)
#    pan4.set_xlabel(r'$q/k_L$')
#    pan4.set_ylabel('synthetic site with maximum population')
    
#    pan4=fig3.add_subplot(2,1,2)    
#    pan4.set_title('Second band')    
#    for i in range(s):    
#        pan4.scatter(kList,[i for j in range(kList.size)],c=pops[:,1,i],vmin=0,vmax=1.0, cmap='Blues', marker='_',linewidths=10)
#    pan4.set_ylabel('Synthetic lattice site')
#    pan5=fig3.add_subplot(3,1,3)    
#    pan5.set_title('Average of First and Second band')    
#    for i in range(s):    
#        pan5.scatter(kList,[i for j in range(kList.size)],c=(pops[:,1,i]+pops[:,0,i])/2.0,vmin=0,vmax=1.0,cmap='Blues', marker='_',linewidths=10)
#    pan5.set_ylabel('Synthetic lattice site')
#    pan5.set_xlabel(r'$k/k_L$')
    
#    fig4=plt.figure()
#    pan4=fig4.add_subplot(1,1,1)
#    pan4.plot(kList,m[:,0],'k-',lw=8)
#    if magCell:
#        for j in range(s-2):
#            pan4.plot(kList,magInUnitCell[:,j])
#    pan4.set_xlabel(r'Crystal momentum [$k_L$]')
#    pan4.set_ylabel('Magnetization')
    if plot:
        cDict={}
        cDict[-2]='m-'
        cDict[-1]='r-'
        cDict[0]='g-'
        cDict[1]='b-'
        cDict[2]='c-'
    #
        fig5=plt.figure()
        pan5=fig5.add_subplot(1,1,1)
        for mF in np.arange(-S,S+1):
            pan5.plot(kList/2.0,pops[:,0,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
        pan5.set_xlabel(r'$q_x$ [$2k_L$]')
        pan5.set_xlim([-0.5,0.5])
        pan5.set_ylabel('Fractional populations')
        pan5.set_title('Lowest band')
        plt.legend()
    
        fig6=plt.figure()
        pan6=fig6.add_subplot(1,1,1)
        for mF in np.arange(-S,S+1):
            pan6.plot(kList/2.0,pops[:,1,mF+S],cDict[mF],label=r'$m_F$='+str(mF))
        pan6.set_xlabel(r'$q_x$ [$2k_L$]')
        pan6.set_xlim([-0.5,0.5])
        pan6.set_ylabel('Fractional populations')
        pan6.set_title('Second band')
        plt.legend()    
#    
    if save:
        filename='SynDimBandStructure_F'+str(S)+'_n'+str(n)+'_Chern'+str(int(c))
        if Flat:
            filename=filename+'Flat'
        print filename
        np.savez(filename, omega=omega,delta=delta,epsilon=epsilon,U=U,kList=kList,E=E,m=m,pops=pops,popAvg=popAvg)
#    
#    a=np.arange(2*S+1)
#    nbar=np.dot(m1pops,a)
#    fig7=plt.figure()
#    pan7=fig7.add_subplot(1,1,1)
#    pan7.plot(kList,nbar)
#    pan7.set_xlabel(r'Crystal momentum [$k_L$]')
#    pan7.set_ylabel('Center of mass')
#    plt.legend()   
    return kList,E,pops[:,0,:], m
    
def plotSynDimBandStructGenAllE(omega, delta, epsilon, U, n,S,m0,kList=np.linspace(-1.0,1.0,600)):
    i=0  
    s=2*S+1
    N=2*n+1
    E=np.zeros((kList.size,s*N))
    m=np.zeros((kList.size,s*N))
    pops=np.zeros((kList.size,s*N,s))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
    for k in kList:
        H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(N*s,N,s)
        popsB=np.sum(evA*np.conjugate(evA),axis=1)
        pops[i]=popsB
        
        m[i]=np.dot(popsB,mFvect)
        i=i+1
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
    for i in range(s):   
        d=panel.scatter(kList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='_')
    panel.set_xlabel(r'$q/k_L$')
    panel.set_ylabel(r'$E/E_L$')
    plt.colorbar(d)
    
    fig3=plt.figure()
    fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
    pan3=fig3.add_subplot(1,1,1)     
    for i in range(s):    
        pan3.scatter(kList,[i for j in range(kList.size)],c=pops[:,0,i],vmin=0,vmax=1.0, cmap='Blues', marker='_',linewidths=10)
    
    pan3.set_ylabel('Synthetic lattice site')
    pan3.set_xlabel(r'Crystal momentum [$k_L$]')
    
def plotPopsOfOmega(omegaMin,omegaMax,step,delta,epsilon,U,n,S,m0,k=0.0):
    i=0  
    s=2*S+1
    N=2*n+1
    omegaList=np.arange(omegaMin,omegaMax,step)
    E=np.zeros((omegaList.size,s*N))
    m=np.zeros((omegaList.size,s*N))
    pops=np.zeros((omegaList.size,s*N,s))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
    for omega in omegaList:
        H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(N*s,N,s)
        popsB=np.sum(evA*np.conjugate(evA),axis=1)
        pops[i]=popsB
        
        m[i]=np.dot(popsB,mFvect)
        i=i+1
    
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
    fig5=plt.figure()
    fig5.suptitle(r'F= '+str(S)+r', $k$ = '+str(np.round(k,2))+r'$k_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
    pan5=fig5.add_subplot(1,1,1)
    for mF in np.arange(-S,S+1):
        pan5.plot(omegaList,pops[:,0,mF+S],cDict[mF], label=r'$m_F$='+str(mF))
    pan5.set_xlabel(r'$\Omega$ [$E_L$]')
    pan5.set_ylabel('Fractional populations')
    #pan5.set_title('Lowest band')
    plt.legend()
    
def plotBandSuperpositionOfT(tlist,omega,delta,epsilon,U,n,S,m0,k=0.0,f2=0.1,f3=0.0,f4=0.0,f5=0.0,plot=True):
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
 #   tlist=np.linspace(0.0,tmax*Erecoil/hbar,num=100)
    
    H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)

    Energy, V = LA.eig(H)
    sort = np.argsort(Energy)
    Esorted, eVsorted = Energy[sort], V[:,sort]
    psi0=np.sqrt(1-f2-f3)*eVsorted[:,0]+np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]
    if S==2:
        psi0=np.sqrt(1-f2-f3-f4-f5)*eVsorted[:,0]+np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]+np.sqrt(f4)*eVsorted[:,3]+np.sqrt(f5)*eVsorted[:,4]
    
    mat=np.array([np.identity(Energy.size)]*tlist.size)

    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    V=np.array([V]*tlist.size)
    Vinv=np.array([Vinv]*tlist.size)
    psi0=np.array([psi0]*tlist.size)


    Uprop = np.einsum('ij,ijk->ijk',np.exp(-1j*np.outer(tlist, Energy)),mat)

    a = np.einsum('ijk,ik->ij',Vinv,psi0)
    b = np.einsum('ijk,ik->ij',Uprop,a)                                                           
    psi = np.einsum('ijk,ik->ij',V,b)
    
    pops=np.absolute(psi)**2.0   
    spinPops=np.sum(pops.reshape(tlist.size,2*n+1,2*S+1),axis=1).transpose()
    
    if plot:
            
        figure=plt.figure()
        panel=figure.add_subplot(1,1,1)
        panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$, f2,f3,f4,f5 = %.2f,%.2f,%.2f,%.2f' %(omega,delta,U,f2,f3,f4,f5))
        for s in range(2*S+1):
            panel.plot(tlist*hbar*1.0e3/Erecoil,spinPops[s],cDict[int(s-S)],label=r'$m_F$='+str(int(s-S)))
    
        panel.set_xlabel('Hold time [ms]')
        plt.legend()
    return spinPops.transpose().flatten()

def plotBandSuperpositionOfTandD(tlist,omega,deltaList,epsilon,U,n,S,m0,phi,k=0.0,f2=0.1,f3=0.0,f4=0.0,f5=0.0,plot=True):
    tlist=np.sort(tlist)    
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
 #   tlist=np.linspace(0.0,tmax*Erecoil/hbar,num=100)
    
    H,kstates=RamanLatHam(k, omega, deltaList[0], epsilon, U, n,S,m0)

    Energy, V = sLA.eigh(H)
    sort = np.argsort(Energy)
    Esorted, eVsorted = Energy[sort], V[:,sort]
    psi0=np.sqrt(1-f2-f3)*eVsorted[:,0]+np.exp(1j*phi)*np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]
    if S==2:
        psi0=np.sqrt(1-f2-f3-f4-f5)*eVsorted[:,0]+np.exp(1j*phi)*np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]+np.sqrt(f4)*eVsorted[:,3]+np.sqrt(f5)*eVsorted[:,4]
    
    psiList=np.zeros((tlist.size,psi0.size),dtype=complex)
    for ind,t in enumerate(tlist):
        if ind==0:
            tpast=0
        else:
            tpast=tlist[ind-1]
        psiList[ind]=psi0
        H,kstates=RamanLatHam(k, omega, deltaList[ind], epsilon, U, n,S,m0)
        Energy, V = sLA.eigh(H)
        V = V + 1j*0.0
        Vinv = np.conjugate(np.transpose(V))
        Uprop = np.exp(-1j*(t-tpast)*Energy)
        a=np.dot(Vinv,psi0)
        b=Uprop*a
        psi0=np.dot(V,b)
        
    
    pops=np.absolute(psiList)**2.0   
    spinPops=np.sum(pops.reshape(tlist.size,2*n+1,2*S+1),axis=1).transpose()
    
    if plot:
            
        figure=plt.figure()
        panel=figure.add_subplot(1,1,1)
        panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$, f2,f3,f4,f5 = %.2f,%.2f,%.2f,%.2f' %(omega,delta,U,f2,f3,f4,f5))
        for s in range(2*S+1):
            panel.plot(tlist*hbar*1.0e3/Erecoil,spinPops[s],cDict[int(s-S)],label=r'$m_F$='+str(int(s-S)))
    
        panel.set_xlabel('Hold time [ms]')
        plt.legend()
    return 
    
def plotBandSuperpositionOfTandDforFit(tlist,phi1,f2,A,phi2,offset):
    tlist=np.sort(tlist)    
    deltaList=sine60(tlist*hbar/Erecoil,A,phi2,offset)
    cDict={}
    cDict[-2]='m-'
    cDict[-1]='r-'
    cDict[0]='g-'
    cDict[1]='b-'
    cDict[2]='c-'
 #   tlist=np.linspace(0.0,tmax*Erecoil/hbar,num=100)
    
    H,kstates=RamanLatHam(k, omega, deltaList[0], epsilon, U, n,S,m0)

    Energy, V = sLA.eigh(H)
    sort = np.argsort(Energy)
    Esorted, eVsorted = Energy[sort], V[:,sort]
    psi0=np.sqrt(1-f2-f3)*eVsorted[:,0]+np.exp(1j*phi1)*np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]
    if S==2:
        psi0=np.sqrt(1-f2-f3-f4-f5)*eVsorted[:,0]+np.exp(1j*phi1)*np.sqrt(f2)*eVsorted[:,1]+np.sqrt(f3)*eVsorted[:,2]+np.sqrt(f4)*eVsorted[:,3]+np.sqrt(f5)*eVsorted[:,4]
    
    psiList=np.zeros((tlist.size,psi0.size),dtype=complex)
    for ind,t in enumerate(tlist):
        if ind==0:
            tpast=0
        else:
            tpast=tlist[ind-1]
        psiList[ind]=psi0
        H,kstates=RamanLatHam(k, omega, deltaList[ind], epsilon, U, n,S,m0)
        Energy, V = sLA.eigh(H)
        V = V + 1j*0.0
        Vinv = np.conjugate(np.transpose(V))
        Uprop = np.exp(-1j*(t-tpast)*Energy)
        a=np.dot(Vinv,psi0)
        b=Uprop*a
        psi0=np.dot(V,b)
        
    
    pops=np.absolute(psiList)**2.0   
    spinPops=np.sum(pops.reshape(tlist.size,2*n+1,2*S+1),axis=1).transpose()
    
    return spinPops.transpose().flatten()
#omega=0.5
#delta=0.0
#epsilon=0.02
#U=4.4
#n=3
#S=2
#m0=0
#
#f3=0.0
#f4=0.0
#f5=0.0
#dataFile=np.load('22Mar2017_noKick.npz')
##
##  
##signalGood=dataFile['signalGood']
##imbal=dataFile['imbalArray']
##cutoff=0.25
##fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
#fieldGoodArray=np.ones(dataFile['tlist'].size,dtype=bool)
#toffs=0.0
#tList=dataFile['tlist'][fieldGoodArray]+toffs
#sort=np.argsort(tList)
#tList=tList[sort]
##
#fractionP=dataFile['fractionP'][fieldGoodArray][sort]
#fraction0=dataFile['fraction0'][fieldGoodArray][sort]
#fractionM=dataFile['fractionM'][fieldGoodArray][sort]
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
##tList=waveDict['pulseDelay']#np.load('SynDimCalibData31Mar16.npz')['tList']
#tRecoils = np.array(tList*Erecoil/hbar)
#fractions=np.array(fractions)
#frac0=np.array(fractions[:,0])
##print tRecoils
##psi0=np.array([0+1j*0.0,1.0+1j*0.0,0.0+1j*0.0])
#rampOnt=0.0000003*Erecoil/hbar
##
##def constrainedBandSup( omega,delta,epsilon,U, n,S,m0):
##    return lambda t,f2: np.array(plotBandSuperpositionOfT(t,omega,delta,epsilon,U,n,S,m0,f2=f2,plot=False))
##BandSupofT=constrainedBandSup(omega,delta,epsilon,U, n,S,m0)
####print H(1.0,4.0,0.0,0.0)
###a=np.array(tRecoils)
###b=np.array(frac0)
##
#phi1guess=0.3*np.pi
#f2guess=0.01
#Aguess=-0.2
#phi2guess=1.2
#offsetguess=0.176
#(phi1,f2,A,phi2,offset),pcov = optimize.curve_fit(plotBandSuperpositionOfTandDforFit,tRecoils[:maxInd],fractions.flatten(), p0=(phi1guess,f2guess,Aguess,phi2guess,offsetguess))
#print (phi1,f2,A,phi2,offset), np.sqrt(np.diag(pcov))
##popt=(0.5,-0.2)
#tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),100)
#pops_fitted=plotBandSuperpositionOfTandDforFit(tForFit,phi1,f2,A,phi2,offset)
##sort=np.argsort(tRecoils)
##tSorted=tRecoils[sort]
#
#pops_fitted=pops_fitted.reshape(tForFit.size,2*S+1).transpose()
#
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
##panel.set_title( r'$\Omega$ = ' + str(np.round(omega,2)) + r' $E_L/h$,$\delta$ = '+str(np.round(delta,3))+r' $E_L/h$, f2 = '+str(np.round(popt[0],3)))#, U= '+str(np.round(popt[2],3))+r'$E_L$')
#
#panel.plot(tList*1.0e3,fractionP,'bo', label=r'$m_F$=+1')
#panel.plot(tList*1.0e3,fraction0,'go', label=r'$m_F$=0')
#panel.plot(tList*1.0e3,fractionM,'ro', label=r'$m_F$=-1')
#
#
#panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[S],'g-')
#panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[S+1],'b-')
#panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[S-1],'r-')
#
#if S==2:
#    panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[S+2],'c-')
#    panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[S-2],'m-')
#    panel.plot(tList*1.0e3,fractionM2,'mo', label=r'$m_F$=-2')
#    panel.plot(tList*1.0e3,fractionP2,'co', label=r'$m_F$=+2')
#panel.set_xlabel('Lattice pulse time [ms]')
#legend()
