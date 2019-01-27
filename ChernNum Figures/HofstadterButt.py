# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 17:35:04 2019

@author: swooty
"""
import numpy as np
from scipy import linalg as sLA
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=3.0
c=4.0/q#0.0#1064.0/790.0#


def TBham2dV2(tx,tm,p,q,kj,km):
    kj=np.array(kj)
    km=np.array(km)
#    H=np.zeros((kj.size,km.size,q,q), dtype=complex)
#    kList = np.zeros((kj.size,q))
    cent=np.int(q/2)
    kList = (np.ones((q,kj.size))*kj).transpose() - 2.0*np.pi*p*(np.arange(q)*np.ones((kj.size,q))-cent)/q
    
    H = np.einsum('ij,ljm->iljm',-2.0*tx*np.cos(kList),np.array([np.diag(np.ones(q))]*km.size))
#    H = np.diag(-2.0*tx*np.cos(kList)) 
    if q>2:
        H = H - tm*np.array([[np.diag(np.ones(q-1),1)]*km.size]*kj.size)
        H = H - tm*np.array([[np.diag(np.ones(q-1),-1)]*km.size]*kj.size)
    try:
        H = H - np.einsum('il,jk->iljk',np.array([tm*np.exp(-1.0j*km*q)]*kj.size),np.diag(np.array([1.0]),-(q-1)))
        H = H - np.einsum('il,jk->iljk',np.array([tm*np.exp(1.0j*km*q)]*kj.size),np.diag(np.array([1.0]),(q-1)))
    except ValueError:
        print 'going to exception'
        H = H - tm*np.exp(-1.0j*km*q)*np.diag(np.array([1.0]),-(q-1))
        H = H - tm*np.exp(1.0j*km*q)*np.diag(np.array([1.0]),(q-1))
        
#    H = H - np.diag(tm*np.exp(1.0j*km*q),q) 
#    
#    for i in np.arange(q):
#        kList[i]=kj-2.0*np.pi*p*(i-cent)/q
#        H[i,i]=-2.0*tx*np.cos(kList[i])
#        if i==0:
#            if q==1:
#                H[i,i]+=-2.0*tm*np.cos(km*q)
#            else:
#                H[i,i+1]=-tm
#                H[i,q-1]=-tm*np.exp(1.0j*km*q)
#        elif i==q-1:
#            H[i,0]=-tm*np.exp(-1.0j*km*q)
#            H[i,i-1]=-tm
#        else:
#            H[i,i+1]=-tm
#            H[i,i-1]=-tm

    return H
    
def getEigenspectrum(tx,tm,p,q, plot=True):
    kjList=np.linspace(-np.pi,np.pi,300)
    kmList=np.linspace(-np.pi/q,np.pi/q,300)
    Egrid=np.zeros((kjList.size,kmList.size,q))
    Vmag=np.zeros((kjList.size,kmList.size,q))
    mgrid=np.zeros((kjList.size,kmList.size))
    cent=np.int(q/2)
    mList=np.arange(q)-cent
    for kjind,kj in enumerate(kjList):
        for kmind, km in enumerate(kmList):
            H,kList=TBham2d(tx,tm,p,q,kj,km)
            E,V = sLA.eigh(H)
            Egrid[kjind,kmind,:] = E
            Vmag[kjind,kmind,:] = V[:,0]*np.conjugate(V[:,0])
            mgrid[kjind,kmind] = np.dot(mList,Vmag[kjind,kmind])

    if plot:        
        figure=plt.figure()
        pan=figure.add_subplot(111)
        pan.imshow(Egrid[:,:,0], cmap='Greys', extent=(kmList[0],kmList[-1],kjList[0],kjList[-1]))
    return Egrid,kjList,kmList,mgrid, Vmag

def getEigenspectrumAll(tx,tm,p,q):
    kjList=np.linspace(-np.pi,np.pi,30)
    kmList=np.linspace(-np.pi/q,np.pi/q,30)
    H = TBham2dV2(tx,tm,p,q,kjList,kmList)

    E,V = np.linalg.eigh(H)
            
    return E

def gcd_recursive(a, b):
    if b == 0:
        return a
    else:
        return gcd_recursive(b, a % b)
    
def mutuallyPrime(q,p):
    return gcd_recursive(p,q)==1
    
qList = np.arange(1,35,1)
fluxList = np.zeros((0))
ElistAll = np.zeros((0))
for qind, q in enumerate(qList):
    pList = np.arange(q-1)+1
    pList = np.append(1,pList[[mutuallyPrime(p,q) for p in pList]])
    for pind, p in enumerate(pList):
        flux = np.float(p)/np.float(q)
        Egrid = getEigenspectrumAll(0.1,0.1,p,q)
        Elist = Egrid.flatten()
        fluxList = np.append(fluxList,[flux]*Elist.size)
        ElistAll = np.append(ElistAll,Elist)
        
np.savez('HofstadtersButt',fluxList=fluxList,ElistAll=ElistAll)
        
fig = plt.figure()
pan = fig.add_subplot(111)
pan.plot(fluxList,ElistAll,'k.', markersize=1)