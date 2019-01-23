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
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=3.0
c=-1064.0/790.0#4.0/q#0.0#
periodic=False
Flat=True
saveBandStruct=False

#
#S=4
#omegaList=np.linspace(0.0001,5.0)
#chernList=np.zeros(omegaList.size)
#chernList2=np.zeros(omegaList.size)
#for ind, omega in enumerate(omegaList):
#    omega2=omega
#    (chernList[ind],chernList2[ind])=plotSynDimBandStructGenTB(omega, 0.0, 0.0, 0.35, S, 0,plots=False,kList=np.linspace(-1.0,1.0,1200))
#    
#fig=plt.figure()
#pan=fig.add_subplot(111)
#pan.plot(omegaList,chernList,label='lowest band')
#pan.plot(omegaList,chernList2,label='lowest bands')
#pan.set_xlabel(r'Raman coupling strength [$E_L$]')
#pan.set_ylabel(r'Chern number from slope')
#pan.set_title(r'F=2, $t_x$= 0.35 $E_L$, $\delta$=$\epsilon$=0.0, cyclic coupling')
    
    
    
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
                F[i,j]=1.0/np.sqrt(2.0)
    return F

def RamanLatHamTB(k, omega, delta, epsilon, tx, S, m0, omega2=0.0,c=c,Flat=True):

    Nspin=int(2*S+1)
    Ntot=Nspin
    Kinetic=np.zeros(Ntot)
    kStates=np.zeros(Ntot)

    for i in range(Ntot):
        spinI=np.float(i-S)
        Kinetic[i]=-tx*np.cos((k+2.0*c*(spinI-m0))*(np.pi))
    H=np.diag(Kinetic)
    H+=delta*Fz(S)
    try:
        H+=((-1.0)**(S+1))*epsilon*(Fz(S)**2.0)
    except ValueError:
        H+=epsilon*(Fz(S)**2.0)
    if Flat:
        H+=(np.sqrt(2.0)/2.0)*omega*FxFlat(S)  
    else:
         H+=(np.sqrt(2.0)/2.0)*omega*Fx(S)
         
    if periodic:
        H[0,Ntot-1]=(np.sqrt(2.0)/2.0)*omega2*1.0/np.sqrt(2.0)
        H[Ntot-1,0]=(np.sqrt(2.0)/2.0)*omega2*1.0/np.sqrt(2.0)
    return H, kStates
    
def plotSynDimBandStructGenTB(omega, delta, epsilon, tx, S, m0,c=c,kList=np.linspace(-1.0,1.0,600),save=False,plots=True):
    i=0  
        
    s=int(2*S+1)
    E=np.zeros((kList.size,s))
    m=np.zeros((kList.size,s))
    pops=np.zeros((kList.size,s,s))
    ktot=np.zeros((kList.size,s))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
  #  mFvect=np.arange(s)
    magInUnitCell=np.zeros((kList.size,s-2))
    for k in kList:
        H,kstates=RamanLatHamTB(k, omega, delta, epsilon, tx, S,m0,c=c)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(s,s)
        popsA=evA*np.conjugate(evA)
        pops[i]=popsA
        
        evB=ev.reshape(s,s)
        popsB=evB*np.conjugate(evB)
        ktot[i]=np.einsum('i,ji->j',kstates,popsB)
        
        m[i]=np.dot(popsA,mFvect)
        for j in range(s-2):
            magInUnitCell[i,j]=(-popsA[0,j]+popsA[0,j+2])/(popsA[0,j]+popsA[0,j+1]+popsA[0,j+2])
        i=i+1
   
    if plots:
        figure=plt.figure()
        panel=figure.add_subplot(1,1,1)
      #  panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        for i in range(s):   
            d=panel.scatter(kList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='_')
        cbar = figure.colorbar(d,ticks=np.array([-S,0,S]))
      #  cbar.ax.tick_params(labelsize=15) 
        cbar.set_label(r'$\langle m \rangle$', size=16)
        panel.set_xlabel(r'Crystal momentum [$k_L$]',size=16)
        panel.set_ylabel(r'Energy [$E_L$]',size=16)
        if saveBandStruct:
            plt.savefig('Z:/My Documents/papers etc/talks and posters/largeBS.jpg',dpi=400)

    
#    fig2=plt.figure()
#    pan2=fig2.add_subplot(1,1,1)    
#    pan2.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
#
#    for i in range(s):   
#        pan2.plot(kList,ktot[:,i],label='band '+str(i))
#    pan2.set_xlabel(r'$q/k_L$')
#    pan2.set_ylabel(r'Total momentum $[k_L]$')
#    plt.legend()
    
    
    kM=kList[np.where(pops[:,0,S-1]==np.max(pops[:,0,S-1]))]
    k0=kList[np.where(pops[:,0,S]==np.max(pops[:,0,S]))]
    kP=kList[np.where(pops[:,0,S+1]==np.max(pops[:,0,S+1]))]
    print kM, k0, kP
    
    chern=2.0*2.0/(3.0*(kP-kM))
    print 'chern number from slope lowest band = ' +str(chern)
    mM=m[kList.size/6,0]
    mP=m[-kList.size/6-1,0]
    print mM, mP
    print 'chern number from displacement lowest band = ' + str((mP-mM)/2.0)
    
    if plots:
        fig3=plt.figure()
    #    fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest band')    
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=pops[:,0,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
        
        pan3.set_ylabel('Synthetic lattice site')
        pan3.set_xlabel(r'Crystal momentum [$k_L$]')
        
        fig3=plt.figure()
#        fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest bands')  
    
    popAvg=np.zeros(pops[:,0,:].shape)
    for ind,k in enumerate(kList):
        for i in range(s):
            popAvg[ind,i]=np.average(pops[ind,0:int(np.ceil(s/q)),i])
    if plots:   
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=popAvg[:,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
    
    kM=kList[np.where(popAvg[:,S-1]==np.max(popAvg[:,S-1]))]
    k0=kList[np.where(popAvg[:,S]==np.max(popAvg[:,S]))]
    kP=kList[np.where(popAvg[:,S+1]==np.max(popAvg[:,S+1]))]
    #print kM, k0, kP
    
    chern2=2.0*2.0/(3.0*(kP-kM))
    print 'chern number from slope lowest '+str(int(np.ceil(s/q))) +' bands = ' +str(chern2)

    if plots:
        pan3.set_ylabel('Synthetic lattice site')
        pan3.set_xlabel(r'Crystal momentum [$k_L$]')
    
        fig3=plt.figure()
 #       fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest bands no edges')      
    
    popAvg2=np.zeros(pops[:,0,:].shape)
    for ind,k in enumerate(kList):
        for i in range(s):
            popAvg2[ind,i]=np.average(pops[ind,0:int(np.ceil(s/q))-1,i])
            
    if plots:
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=popAvg2[:,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
            
        fig3=plt.figure()
#        fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Edge band')      
            
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=pops[:,int(np.ceil(s/q))-1,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)

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
    mBands=np.dot(popAvg,mFvect)
    
    mM=mBands[kList.size/6]
    mP=mBands[-kList.size/6-1]
    print mM, mP
    print 'chern number from displacement lowest band = ' + str((mP-mM)/2.0)
    
    if plots:
        fig4=plt.figure()
        pan4=fig4.add_subplot(1,1,1)
        pan4.plot(kList,m[:,0],'k-',lw=2, label = 'lowest band')
        pan4.plot(kList,mBands,'b-',lw=2, label = 'lowest '+str(int(np.ceil(s/q))) +' bands' )
    #    for j in range(s-2):
    #        pan4.plot(kList,magInUnitCell[:,j])
        pan4.set_xlabel(r'Crystal momentum [$k_L$]')
        pan4.set_ylabel('Magnetization')
        legend()
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
            pan5.plot(kList,pops[:,0,mF+S], cDict[mF], label=r'$m_F$='+str(mF))
            pan5.plot(kList+2.0,pops[:,1,mF+S],cDict[mF])
        pan5.set_xlabel(r'Crystal momentum [$k_L$]')
        pan5.set_ylabel('Fractional populations')
        pan5.set_title('Lowest band to second band')
        plt.legend()
    #
    #    fig6=plt.figure()
    #    pan6=fig6.add_subplot(1,1,1)
    #    for mF in np.arange(-S,S+1):
    #        pan6.plot(kList,pops[:,1,mF+S],label=r'$m_F$='+str(mF))
    #    pan6.set_xlabel(r'Crystal momentum [$k_L$]')
    #    pan6.set_ylabel('Fractional populations')
    #    pan6.set_title('Second band')
    #    plt.legend()    
    
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
    return kList, E, m, chern, chern2

def plotSynDimBandStructGenTBClean(omega, delta, epsilon, tx, S, m0,c=c,kList=np.linspace(-1.0,1.0,600),Flat=True,save=False,plots=True):
    i=0  
        
    s=int(2*S+1)
    E=np.zeros((kList.size,s))
    m=np.zeros((kList.size,s))
    pops=np.zeros((kList.size,s,s))
    ktot=np.zeros((kList.size,s))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
  #  mFvect=np.arange(s)
    magInUnitCell=np.zeros((kList.size,s-2))
    for k in kList:
        H,kstates=RamanLatHamTB(k, omega, delta, epsilon, tx, S,m0,c=c,Flat=Flat)
        Energies, eigenstates = LA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(s,s)
        popsA=evA*np.conjugate(evA)
        pops[i]=popsA
        
        evB=ev.reshape(s,s)
        popsB=evB*np.conjugate(evB)
        ktot[i]=np.einsum('i,ji->j',kstates,popsB)
        
        m[i]=np.dot(popsA,mFvect)
        for j in range(s-2):
            magInUnitCell[i,j]=(-popsA[0,j]+popsA[0,j+2])/(popsA[0,j]+popsA[0,j+1]+popsA[0,j+2])
        i=i+1
   
    if plots:
        figure=plt.figure()
        panel=figure.add_subplot(1,1,1)
      #  panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        for i in range(s):   
            d=panel.scatter(kList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='_')
        cbar = figure.colorbar(d,ticks=np.array([-S,0,S]))
      #  cbar.ax.tick_params(labelsize=15) 
        cbar.set_label(r'$\langle m \rangle$', size=16)
        panel.set_xlabel(r'Crystal momentum [$k_L$]',size=16)
        panel.set_ylabel(r'Energy [$E_L$]',size=16)
        if saveBandStruct:
            plt.savefig('Z:/My Documents/papers etc/talks and posters/largeBS.jpg',dpi=400)


    
    if plots:
        fig3=plt.figure()
    #    fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest band')    
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=pops[:,0,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
        
        pan3.set_ylabel('Synthetic lattice site')
        pan3.set_xlabel(r'Crystal momentum [$k_L$]')
        
        fig3=plt.figure()
#        fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest bands')  
    
    popAvg=np.zeros(pops[:,0,:].shape)
    for ind,k in enumerate(kList):
        for i in range(s):
            popAvg[ind,i]=np.average(pops[ind,0:int(np.ceil(s/q)),i])
    if plots:   
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=popAvg[:,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
    

    if plots:
        pan3.set_ylabel('Synthetic lattice site')
        pan3.set_xlabel(r'Crystal momentum [$k_L$]')
    
        fig3=plt.figure()
 #       fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Lowest bands no edges')      
    
    popAvg2=np.zeros(pops[:,0,:].shape)
    for ind,k in enumerate(kList):
        for i in range(s):
            popAvg2[ind,i]=np.average(pops[ind,0:int(np.ceil(s/q))-1,i])
            
    if plots:
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=popAvg2[:,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)
            
        fig3=plt.figure()
#        fig3.suptitle(r'F= '+str(S)+r', $\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
        pan3=fig3.add_subplot(1,1,1)    
        pan3.set_title('Edge band')      
            
        for i in range(s):    
            pan3.scatter(kList,[i for j in range(kList.size)],c=pops[:,int(np.ceil(s/q))-1,i],cmap='Blues',vmin=0.0,vmax=1.0/S,marker='_',linewidths=10)


    return kList, E, m, pops