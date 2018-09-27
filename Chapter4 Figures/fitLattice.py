# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 11:55:14 2016

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

def LatHam(U, n, k):
    
    Nlat=2*n+1

    Kinetic=np.zeros(Nlat)

    for i in range(Nlat):
        latI=i-n
        Kinetic[i]=(k-2.0*latI)**2.0
    H=np.diag(Kinetic)

    for i in range(Nlat):
        for j in range(Nlat):
            if np.abs(i-j)==1:
                H[i,j]=-U/4.0         
    return H
    
    
def propagateLatHam(t, U,k=0,psi0='norm' ): 
    n=7
    t=np.array(t)
    H = LatHam(U,n,k)
    Energy, V = LA.eig(H)
    if psi0=='norm':
        psi0 = np.zeros(Energy.size)
        psi0[n]=1.0


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
    pops=pops.reshape(t.size,2*n+1)
    popsThree=pops[:,n-1:n+2]    
    popsThree=popsThree.flatten()             
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    return popsThree, psi
    
def propagateLatHamFit(t, U,k=0,psi0='norm' ): 
    n=7
    t=np.array(t)
    H = LatHam(U,n,k)
    Energy, V = LA.eig(H)
    if psi0=='norm':
        psi0 = np.zeros(Energy.size)
        psi0[n]=1.0


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
    pops=pops.reshape(t.size,2*n+1)
    popsThree=pops[:,n-1:n+2]    
    popsThree=popsThree.flatten()             
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    return popsThree   
    
def plotPulsedPops(tf, U,k=0):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    pops,psi=propagateLatHam(tlist,U,k=k)
    pops=pops.reshape(tlist.size,3)

    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'U= %.2f $E_L$, q= %.2f $k_L$' %(U,k))
    panel.plot(tlist*hbar*1.0e6/Erecoil,pops[:,0],'b-',label=r'$k$='+str(k+2)+r' $k_L$')
    panel.plot(tlist*hbar*1.0e6/Erecoil,pops[:,1],'g-',label=r'$k$='+str(k)+r' $k_L$' )
    panel.plot(tlist*hbar*1.0e6/Erecoil,pops[:,2],'r-',label=r'$k$='+str(k-2)+r' $k_L$')
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    return

def PulseTrainPops(Npulses,pulseTime, U,k=0):
    tRecoils = pulseTime*Erecoil/hbar
    pops,psi=propagateLatHam(tRecoils,U,k=k)
    
    for i in range(Npulses-1):
        pops,psiNew=propagateLatHam(tRecoils,0,k=k,psi0=psi)
        pops,psi=propagateLatHam(tRecoils,U,k=k,psi0=psiNew)

    return pops
    
def plotPulseTrainOfT(tlist,Npulses,U,k=0):
    popsAll=np.zeros((tlist.size,3))
    
    for ind,t in enumerate(tlist):
        popsAll[ind] = PulseTrainPops(Npulses,t,U,k=k)
        
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title('N,U,k = %.0f,%.2f,%.2f' %(Npulses,U,k))
    panel.plot(tlist*1.0e6,popsAll[:,0],'b-',label=r'$k$='+str(k+2)+r' $k_L$')
    panel.plot(tlist*1.0e6,popsAll[:,1],'g-',label=r'$k$='+str(k)+r' $k_L$' )
    panel.plot(tlist*1.0e6,popsAll[:,2],'r-',label=r'$k$='+str(k-2)+r' $k_L$')
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    
    return
    
def plotPulseTrainOfQ(qlist,pulseTime,Npulses,U):
    popsAll=np.zeros((qlist.size,3))
    
    for ind,q in enumerate(qlist):
        popsAll[ind] = PulseTrainPops(Npulses,pulseTime,U,k=q)
        
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title('N,U,pulseTime = %.0f,%.2f,%.6f' %(Npulses,U,pulseTime))
    panel.plot(qlist,popsAll[:,0],'b-',label=r'$k = q + 2k_L$')
    panel.plot(qlist,popsAll[:,1],'g-',label=r'$k = q$')
    panel.plot(qlist,popsAll[:,2],'r-',label=r'$k = q - 2k_L$')
    panel.set_xlabel('Initial momentum [k_L]')
    plt.legend()
    
    return    
    
def plotPulsedPopsOfK(tf, U,kList=np.linspace(-1.0,1.0,100)):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    popsAll=np.zeros((kList.size,tlist.size,3))
    for ind,k in enumerate(kList):
        
        pops, psi=propagateLatHam(tlist,U,k=k)
        pops=pops.reshape(tlist.size,3)
        popsAll[ind]=pops

    fig2=plt.figure()
    fig2.clear()
    gs = gridspec.GridSpec(1,3)
    gs.update(left=0.1, right=0.99, wspace=0.005, hspace=0.0)
    pan=fig2.add_subplot(gs[0])
    pan.imshow(popsAll[:,:,0],extent=(tlist[0]*hbar*1e3/Erecoil,tlist[-1]*hbar*1e3/Erecoil,kList[-1],kList[0]),vmin=0.0,vmax=1.0,aspect='equal')
    pan.set_ylabel(r'Quasimomentum [$k_L$]')
    pan.set_xlabel('Lattice pulse time [ms]')

    pan2=fig2.add_subplot(gs[1])
    pan2.imshow(popsAll[:,:,1],extent=(tlist[0]*hbar*1e3/Erecoil,tlist[-1]*hbar*1e3/Erecoil,kList[-1],kList[0]),vmin=0.0,vmax=1.0,aspect='equal')
   # pan2.set_ylabel(r'Quasimomentum [$k_L$]')
    pan2.set_xlabel('Lattice pulse time [ms]')
        
    pan3=fig2.add_subplot(gs[2])
    a=pan3.imshow(popsAll[:,:,2],extent=(tlist[0]*hbar*1e3/Erecoil,tlist[-1]*hbar*1e3/Erecoil,kList[-1],kList[0]),vmin=0.0,vmax=1.0,aspect='equal')
    #pan3.set_ylabel(r'Quasimomentum [$k_L$]')
    pan3.set_xlabel('Lattice pulse time [ms]')
   # plt.colorbar(a)

    return
    
def adiabatic(U,n,k,plot=True):
    H=LatHam(U,n,k)
    E,Vmat=sLA.eigh(H)
    sort=np.argsort(E)
    V0=Vmat[:,sort[0]]
    xlist=np.linspace(-200.0,200.0,5000)
    klist=np.arange(2*n+1)-n
    basis=np.exp(np.outer(-klist*2.0j,xlist))
    psi=np.dot(basis.transpose(),V0)
    psiMag=psi*np.conjugate(psi)
    wan=psi*np.sin(xlist)/xlist
    
    wanMag=wan*np.conjugate(wan)
    wan=wan/np.sqrt(np.sum(wanMag*(xlist[1]-xlist[0])))
    wanMag=wan*np.conjugate(wan)
    wanFft=np.fft.fft(wanMag)
    freq=np.fft.fftfreq(xlist.size,d=xlist[1]-xlist[0])*2.0*np.pi
    freq=np.fft.fftshift(freq)
    wanFft=np.fft.fftshift(wanFft)
    
    momentumExpectation=np.sum(np.conj(V0)*klist*V0)/np.sum(np.conj(V0)*V0)
    print 'momentuum expectation value = %.3f' %(momentumExpectation)
    if plot:
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.plot(xlist,psiMag,'k-')
        pan.plot(xlist,wanMag,'r-')
        
        fig2=plt.figure()
        pan2=fig2.add_subplot(111)
        pan2.plot(freq,wanFft.real,'ko',freq,wanFft.imag,'bo',freq,wanFft*np.conj(wanFft),'r-')
        
    return V0, psi, xlist
    
def wannier(U,n,plot=True,save=False):
    klist=np.linspace(-1.0,1.0,2000)
    V0,psi,xlist=adiabatic(U,n,0,plot=False)
    psiList=np.zeros((klist.size,xlist.size),dtype=complex)
    for i,k in enumerate(klist):
        V0,psiList[i],xlist=adiabatic(U,n,k,plot=False)
        psiList[i]=np.exp(1.0j*k*xlist)*psiList[i]
    wan=np.sum(psiList,axis=0)*(klist[1]-klist[0])/2.0
    
    wanMag=wan*np.conjugate(wan)
    wanFft=np.fft.fft(wan)
    freq=np.fft.fftfreq(xlist.size,d=(xlist[1]-xlist[0])/(2.0*np.pi))
    freq=np.fft.fftshift(freq)
    wanFft=np.fft.fftshift(wanFft)
    
    if plot:
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.plot(xlist,wanMag,'r-')
        
        fig2=plt.figure()
        pan2=fig2.add_subplot(111)
        pan2.plot(freq,wanFft.real,'ko',freq,wanFft.imag,'bo',freq,wanFft*np.conj(wanFft),'r-')
    if save:
        filename='wannier_U_%.2f' %(U)
        np.savez(filename,U=U,n=n,psiList=psiList,wan=wan,wanMag=wanMag,wanFft=wanFft,wanFftMag=wanFft*np.conj(wanFft))
    return wan, wanFft
    
def blochCollect(U,n,plot=True,save=False):
    qlist=np.linspace(-1.0,1.0,500)
    N=2*n+1
    klist=np.linspace(-N,N,qlist.size*N)
    bC=np.zeros(klist.size)
    for i,q in enumerate(qlist):
        inds=np.arange(N)*qlist.size+i
        H=LatHam(U,n,q)
        E,Vmat=sLA.eigh(H)
        sort=np.argsort(E)
        V0=Vmat[:,sort[0]]
        bC[inds]+=np.flipud(V0)**2.0
    if plot:
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.plot(klist,bC,'k-')
    if save:
        filename='blochCollect_U_%.2f' %(U)
        np.savez(filename,U=U,n=n,klist=klist,bC=bC)
    return 
    
#
#roiList=[[445,475,470,495], [505, 535, 470, 495], [570, 600, 470,495]]
#roiList2=[[450,450,465,465], [500, 550, 455, 505], [500, 550, 505, 555]]
#filestart=34
#filestop=63
#fileroot = 'Y:/Data/2017/October/25/PIXIS_25Oct2017' 
#counts, fractions, waveDict, probeAvg = batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc='right') 
#
#time=waveDict['pulseDelay']
#tRecoils = time*Erecoil/hbar
#fractions=np.array(fractions)
#
#Uguess=2.0
#popt,pcov = optimize.curve_fit(propagateLatHamFit,tRecoils,fractions.flatten(), p0=(Uguess))
#print popt,pcov
#tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),100)
#pops_fitted=propagateLatHamFit(tForFit,*popt)
#pop0 = np.array([pops_fitted[i*3] for i in range(tForFit.size)])
#pop1 = np.array([pops_fitted[i*3+1] for i in range(tForFit.size)]) 
#pop2 = np.array([pops_fitted[i*3+2] for i in range(tForFit.size)]) 
#
#
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.set_title('U = '+str(np.round(popt[0],4))+' E_L')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
#panel.plot(time*1e6,fractions[:,0],'bo', label='kx=+2kR')
#panel.plot(time*1e6,fractions[:,1],'go', label='kx=0')
#panel.plot(time*1e6,fractions[:,2],'ro', label='kx=-2kR')
#panel.plot(tForFit*1e6*hbar/Erecoil,pop0,'b-')
#panel.plot(tForFit*1e6*hbar/Erecoil,pop1,'g-')
#panel.plot(tForFit*1e6*hbar/Erecoil,pop2,'r-')
#panel.set_xlabel('Lattice pulse time [us]')
#legend()
#
#

