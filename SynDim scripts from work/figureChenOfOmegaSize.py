# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:06:41 2017

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import lineFit2
import tightBindingHam as tbh
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

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
    
def parabola(xlist,A,B,x0):
    return A*(xlist-x0)**2.0+B
    
def parabolaLine(xlist,A,B,C,x0):
    return A*(xlist-x0)**2.0+B*xlist+C
    
def gaussian(xlist,sigma,A,x0,B):
    return A*np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))+B
    
def fitParabola(xlist,ylist,yerr,plot=False):
    (A,B,x0), cov=optimize.curve_fit(parabola,xlist,ylist,p0=(-4.0,np.max(ylist),xlist[np.where(ylist==np.max(ylist))[0][0]]),sigma=yerr,absolute_sigma=False)
    (dA,dB,dx0)=np.sqrt(np.diag(cov))
    
#    (A,B,C,x0), cov=optimize.curve_fit(parabolaLine,xlist,ylist,p0=(-4.0,np.max(ylist),0.1,xlist[np.where(ylist==np.max(ylist))[0][0]]),sigma=yerr,absolute_sigma=False)
#    (dA,dB,dC,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        xforplot=np.linspace(np.min(xlist),np.max(xlist),num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.plot(xforplot,parabola(xforplot,A,B,x0))
        pan.set_title('A=%.2f +/-%.3f , B=%.2f +/-%.3f , x0=%.2f +/-%.3f ' %(A,dA,B,dB,x0,dx0))
#        pan.plot(xforplot,parabolaLine(xforplot,A,B,C,x0))
#        pan.set_title('A=%.2f +/-%.3f , B=%.2f +/-%.3f , C=%.2f +/-%.3f, x0=%.2f +/-%.3f ' %(A,dA,B,dB,C,dC,x0,dx0))
        
    return x0,dx0
    
    
def fitGauss(xlist,ylist,yerr,plot=False):
    (sigma,A,x0,B), cov=optimize.curve_fit(gaussian,xlist,ylist,p0=(0.56,np.max(ylist),xlist[np.where(ylist==np.max(ylist))[0][0]],0.0),sigma=yerr,absolute_sigma=False)
    (dsigma,dA,dx0,dB)=np.sqrt(np.diag(cov))
    
#    (A,B,C,x0), cov=optimize.curve_fit(parabolaLine,xlist,ylist,p0=(-4.0,np.max(ylist),0.1,xlist[np.where(ylist==np.max(ylist))[0][0]]),sigma=yerr,absolute_sigma=False)
#    (dA,dB,dC,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        xforplot=np.linspace(np.min(xlist),np.max(xlist),num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.plot(xforplot,gaussian(xforplot,sigma,A,x0,B))
        pan.set_title('sigma=%.2f +/-%.3f , A=%.2f +/-%.3f , x0=%.2f +/-%.3f, B=%.2f +/-%.3f ' %(sigma,dsigma,A,dA,x0,dx0,B,dB))
#        pan.plot(xforplot,parabolaLine(xforplot,A,B,C,x0))
#        pan.set_title('A=%.2f +/-%.3f , B=%.2f +/-%.3f , C=%.2f +/-%.3f, x0=%.2f +/-%.3f ' %(A,dA,B,dB,C,dC,x0,dx0))
        
    return x0,dx0
    
    
def getExpectationValue(xlist,ylist,yerr,plot=False):
   # dx=xlist[1]-xlist[0]

    tot = np.sum(ylist/(yerr**2.0))

    xTot = np.sum(xlist*ylist/(yerr**2.0))

    xExp = xTot/tot

    yAtXexp = np.interp(xExp,xlist,ylist)
    
    dExp = np.sqrt(np.sum((xlist-xExp)**2.0))/tot
    yAtXexpP = np.interp(xExp+dExp,xlist,ylist)
    yAtXexpM = np.interp(xExp-dExp,xlist,ylist)
    
    if plot:
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(xlist,ylist,yerr=yerr,fmt='bo')
        pan.vlines(np.array([xExp,xExp+dExp,xExp-dExp]),np.array([0.0,0.0,0.0]),np.array([yAtXexp,yAtXexpP,yAtXexpM]))
        
    
    return xExp, dExp
    
    
def getMaxFitTheory(qlist,fracList,sigma,indC='nan',plot=False):
    w=np.int(1.0*qlist.size/4.0)
    sort=np.argsort(qlist)
    if indC == 'nan':
        indC=np.where(fracList[sort]==np.max(fracList[sort]))[0][0]
    indCent=int(qlist.size/2)
    shift=indCent-indC
    qlistShift=np.roll(qlist[sort],shift)
    if shift>0:
        qlistShift[0:shift]=qlistShift[0:shift]-2.0
    elif shift<0:
        qlistShift[shift:qlist.size]=qlistShift[shift:qlist.size]+2.0
    qlistShift=qlistShift[indCent-w+1:indCent+w]

    fracShift=np.roll(fracList[sort],shift)[indCent-w+1:indCent+w]
    sigmaShift=np.roll(sigma[sort],shift)[indCent-w+1:indCent+w]
    
    xP,dxP=fitGauss(qlistShift,fracShift,sigmaShift,plot=plot)
 #   xP,dxP=getExpectationValue(qlistShift,fracShift,sigmaShift,plot=plot)
    if xP>1.0:
        xP=xP-2.0
    if xP<-1.0:
        xP=xP+2.0
        
    return xP,dxP
    
def getChern(omega, delta, epsilon, U, n,S, m0,kList=np.linspace(-1.0,1.0,600), plot=True):
    i=0  
        
    s=int(2*S+1)
    N=2*n+1
    E=np.zeros((kList.size,s*N))
    m=np.zeros((kList.size,s*N))
    pops=np.zeros((kList.size,s*N,s))
    ktot=np.zeros((kList.size,s*N))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])

    for k in kList:
        H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
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

        i=i+1
        
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
            if S <=2:
                pan5.plot(kList,pops[:,0,mF+S], cDict[mF], label=r'$m_F$='+str(mF))
            else:
                pan5.plot(kList,pops[:,0,mF+S])#, cDict[mF], label=r'$m_F$='+str(mF))
        pan5.set_xlabel(r'Crystal momentum [$k_L$]')
        pan5.set_ylabel('Fractional populations')
        pan5.set_title('Lowest band')
        plt.legend()
    
    if c>0:
        indM=kList.size/6
        ind0=kList.size/2
        indP=5*kList.size/6
    elif c<0:
        indM=5*kList.size/6
        ind0=kList.size/2
        indP=kList.size/6
    else:
        indM=kList.size/2
        ind0=kList.size/2
        indP=kList.size/2
    xMth,dxMth=getMaxFitTheory(kList,pops[:,0,S-1],np.ones(pops[:,0].size)*0.001,indC=indM,plot=plot)  #[np.where(pops[:,0,S-1]==np.max(pops[:,0,S-1]))]
    x0th,dx0th=getMaxFitTheory(kList,pops[:,0,S],np.ones(pops[:,0].size)*0.001,indC=ind0,plot=plot)
    xPth,dxPth=getMaxFitTheory(kList,pops[:,0,S+1],np.ones(pops[:,0].size)*0.001,indC=indP,plot=plot)
    


    Ath,Bth,dAth,dBth=lineFit2.lineFit(np.array([xMth,x0th,xPth]),np.array([-1,0,1]),'q','max. m',errorbars=True,yerr=np.array([dxMth,dx0th,dxPth]),plot=plot)
    print xMth, x0th, xPth
    print 'chern theory = '+ str(Ath*2/3) + '+/-'+str(dAth*2/3)
   
    chern =  Ath*2.0/3.0    
    return pops, chern
    
    
#q=3.0
#c=4.0/q#1064.0/790.0#0.0#
#Flat=False
#
#omegaList=np.linspace(0.0,2.0,50)
#chernList1=np.zeros(omegaList.size)
#chernList2=np.zeros(omegaList.size)
#for ind, omega in enumerate(omegaList):
#    pops, chernList1[ind]=getChern(omega, 0.0, 0.02, 4.4, 7,1,0,kList=np.linspace(-1.0,1.0,200), plot=False)
#    pops, chernList2[ind]=getChern(omega, 0.0, 0.02, 4.4, 7,2,0,kList=np.linspace(-1.0,1.0,200), plot=False)
#    
#figure=plt.figure()
#pan=figure.add_subplot(111)
#pan.plot(omegaList,chernList1, 'bo')
#pan.set_xlabel(r'$\Omega$ [$E_L$]')
#pan.set_ylabel(r'Chern number')
#
#figure2=plt.figure()
#pan2=figure2.add_subplot(111)
#pan2.plot(omegaList,chernList2, 'bo')
#pan2.set_xlabel(r'$\Omega$ [$E_L$]')
#pan2.set_ylabel(r'Chern number')
#filename = 'ChernOfOmegaNotFlat'
#np.savez(filename, c=c,Flat=Flat, omegaList = omegaList, chernList1=chernList1,chernList2=chernList2, 
#         delta=0.0,epsilon=0.02,U=4.4,m0=0,kList=np.linspace(-1.0,1.0,200))

#q=3.0
#c=4.0/q#1064.0/790.0#0.0#
#Flat=True
#Slist=np.arange(1,21)
#chernList=np.zeros(Slist.size)
#for ind, Sloc in enumerate(Slist):
#    pops, chernList[ind]=getChern(0.5, 0.0, 0.0, 4.4, 7,Sloc,0,kList=np.linspace(-1.0,1.0,200), plot=False)
#    print Sloc
#    
#    
#figure3=plt.figure()
#pan3=figure3.add_subplot(111)
#pan3.plot(Slist,chernList, 'bo')
#pan3.set_xlabel(r'F')
#pan3.set_ylabel(r'Chern number')
#
#filename = 'ChernOfSystemSizeOmega05'
#np.savez(filename, c=c,Flat=Flat, Slist = Slist, chernList=chernList, omega=0.5,delta=0.0,epsilon=0.0,U=4.4,m0=0,kList=np.linspace(-1.0,1.0,200))