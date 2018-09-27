# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 16:24:12 2017

@author: dng5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import linalg as sLA

date='28Feb2017'
posFile=np.load(date+'_posKick.npz')
negFile=np.load(date+'_negKick.npz')
plotSimul=True
S=2

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
c=1064.0/790.0#0.0#
saveFileName=date+'_F'+str(S)+'_chern_'+str(int(c))
print saveFileName
Flat=False

if S==2:
    weights=np.array([1.006,1.0,1.034,0.992,0.948])
else:
    weights=np.array([1.0,1.045,1.077])
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
    if Flat:
        H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[FxFlat(S)]*Nlat)    
    else:
         H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
    for i in range(Ntot):
        for j in range(Ntot):
            if np.abs(i-j)==(2*S+1):
                H[i,j]=U/4.0         
    return H, kStates
    
#if plotSimul:
omega=0.5 
delta=0.0
epsilon=0.02
U=4.4
n=7
m0=0
s=2*S+1
N=2*n+1
kList=np.linspace(-1.0,1.0,300)

def getPops(omega,delta,epsilon,U,n,m0,kList,getChern=False):

    E=np.zeros((kList.size,s*N))
    m=np.zeros((kList.size,s*N))
    mShort=np.zeros((kList.size,s*N))
    pops=np.zeros((kList.size,s*N,s))
    ktot=np.zeros((kList.size,s*N))

    mFvect=np.array([np.float(j-S) for j in range(s)])
    mFshort=np.zeros(mFvect.shape)
    mFshort[S-1]=-1.0
    mFshort[S+1]=1.0
    for i,k in enumerate(kList):
        H,kstates=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
        Energies, eigenstates = sLA.eig(H)
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
        mShort[i]=np.dot(popsA,mFshort)/np.sum(popsA[S-1:S+2,0])
        i=i+1
        
    if getChern:
        kM=kList[np.where(pops[:,0,S-1]==np.max(pops[:,0,S-1]))]
        kP=kList[np.where(pops[:,0,S+1]==np.max(pops[:,0,S+1]))]
        
        chern=2.0*2.0/(3.0*(kP-kM))
        print 'chern number from simul via slope= ' +str(chern)
        
        m0=(mShort[np.divide(kList.size,2),0]+mShort[np.divide(kList.size,2)-1,0])/2
        mP=(mShort[5*np.divide(kList.size,6),0]+mShort[5*np.divide(kList.size,6)-1,0])/2
        mM=(mShort[np.divide(kList.size,6),0]+mShort[np.divide(kList.size,6)-1,0])/2
        chern2=np.average([m0-mM,mP-m0])
        print 'chern number from simul via magnetization= ' +str(chern2)
    return pops
    
        
        
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

    
fractionP=np.append(posFile['fractionP'],negFile['fractionP'])
fraction0=np.append(posFile['fraction0'],negFile['fraction0'])
fractionM=np.append(posFile['fractionM'],negFile['fractionM'])
weightSum=fractionP*weights[S+1]+fraction0*weights[S]+fractionM*weights[S-1]

if S==2:
    fractionP2=np.append(posFile['fractionP2'],negFile['fractionP2'])
    fractionM2=np.append(posFile['fractionM2'],negFile['fractionM2'])
    weightSum=weightSum+fractionP2*weights[S+2]+fractionM2*weights[S-2]
    fractionP2w=fractionP2*weights[S+2]/weightSum
    fractionM2w=fractionM2*weights[S-2]/weightSum
fractionPw=fractionP*weights[S+1]/weightSum
fraction0w=fraction0*weights[S]/weightSum
fractionMw=fractionM*weights[S-1]/weightSum  
    
qlist=np.append(posFile['qlist'],negFile['qlist'])
qlistFit=np.append(posFile['qlistFit'],negFile['qlistFit'])
inRange=((qlistFit<1.0)&(qlistFit>-1.0))
numPoints=np.append(posFile['numPoints'],negFile['numPoints'])
sigmaP=np.append(posFile['sigmaP'],negFile['sigmaP'])/np.sqrt(numPoints)
sigma0=np.append(posFile['sigma0'],negFile['sigma0'])/np.sqrt(numPoints)
sigmaM=np.append(posFile['sigmaM'],negFile['sigmaM'])/np.sqrt(numPoints)
if S==2:
    sigmaP2=np.append(posFile['sigmaP2'],negFile['sigmaP2'])/np.sqrt(numPoints)
    sigmaM2=np.append(posFile['sigmaM2'],negFile['sigmaM2'])/np.sqrt(numPoints)



fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 
#pan1.errorbar(qlist,fractionP2,yerr=sigmaP2/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
pan1.errorbar(qlistFit[inRange],fractionPw[inRange],yerr=sigmaP[inRange], fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(qlistFit[inRange],fraction0w[inRange],yerr=sigma0[inRange], fmt='go', label=r'$m_F$=0')
pan1.errorbar(qlistFit[inRange],fractionMw[inRange],yerr=sigmaM[inRange], fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan1.errorbar(qlistFit[inRange],fractionP2w[inRange],yerr=sigmaP2[inRange], fmt='co', label=r'$m_F$=+2')
    pan1.errorbar(qlistFit[inRange],fractionM2w[inRange],yerr=sigmaM2[inRange], fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
#if plotSimul:
pops=getPops(omega,delta,epsilon,U,n,m0,kList,getChern=True)
popsM=getPops(omega,-0.01,epsilon,U,n,m0,kList)
popsP=getPops(omega,0.01,epsilon,U,n,m0,kList)
for mF in np.arange(-S,S+1):
    pan1.plot(kList,pops[:,0,mF+S],cDict[mF])
    pan1.fill_between(kList,popsM[:,0,mF+S],popsP[:,0,mF+S],facecolor=cDict[mF][:-1],alpha=0.3)
   # pan1.plot(kList,popsP[:,0,mF+S],cDict[mF]+'.')
pan1.set_xlabel('quasimomentum [k_L]')
pan1.set_ylabel('Fractional populations')
#legend()

fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(qlistFit[inRange],[1 for j in range(qlistFit[inRange].size)],c=fractionPw[inRange],s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(qlistFit[inRange],[0 for j in range(qlistFit[inRange].size)],c=fraction0w[inRange],s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(qlistFit[inRange],[-1 for j in range(qlistFit[inRange].size)],c=fractionMw[inRange],s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
if S==2:
    pan2.scatter(qlistFit[inRange],[2 for j in range(qlistFit[inRange].size)],c=fractionP2w[inRange],s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
    pan2.scatter(qlistFit[inRange],[-2 for j in range(qlistFit[inRange].size)],c=fractionM2w[inRange],s=30,vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)

pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Crystal momentum [$k_L$]')

qM=qlist[np.where(fractionM==np.max(fractionM))]
q0=qlist[np.where(fraction0==np.max(fraction0))]
qP=qlist[np.where(fractionP==np.max(fractionP))]

mag=(fractionP[inRange]-fractionM[inRange])/(fractionP[inRange]+fraction0[inRange]+fractionM[inRange])
mP=mag[np.where(np.abs(qlistFit[inRange]-2.0/3.0)==np.min(np.abs(qlistFit[inRange]-2.0/3.0)))]
print mP
m0=mag[np.where(np.abs(qlistFit[inRange])==np.min(np.abs(qlistFit[inRange])))]
print m0
mM=mag[np.where(np.abs(qlistFit[inRange]+2.0/3.0)==np.min(np.abs(qlistFit[inRange]+2.0/3.0)))]
print mM
chernD2=np.average([m0-mM,mP-m0])
chernD=2.0*2.0/(3.0*(qP-qM))
print 'chern number measured via slope = ' +str(chernD)
print 'chern number measured via magnetization = ' +str(chernD2)    

np.savez(saveFileName,qlist=qlist[inRange],fractionP=fractionP[inRange],fraction0=fraction0[inRange],weights=weights,
         fractionM=fractionM[inRange],sigmaP=sigmaP[inRange],sigma0=sigma0[inRange],
         sigmaM=sigmaM[inRange], qlistSimul=kList,pops=pops[:,0,:2*S+1],omega=omega,delta=delta,epsilon=epsilon,U=U,qlistFit=qlistFit)
         
if S==2:
    np.savez(saveFileName,qlist=qlist[inRange],fractionP2=fractionP2[inRange],fractionP=fractionP[inRange],
             fraction0=fraction0[inRange],fractionM=fractionM[inRange],weights=weights,
             fractionM2=fractionM2[inRange],sigmaP2=sigmaP2[inRange],sigmaP=sigmaP[inRange],sigma0=sigma0[inRange],
             sigmaM=sigmaM[inRange], sigmaM2=sigmaM2[inRange], qlistSimul=kList, 
            pops=pops[:,0,:2*S+1],omega=omega,delta=delta,epsilon=epsilon,U=U,qlistFit=qlistFit)