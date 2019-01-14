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


def TBham2d(tx,tm,p,q,kj,km):

    H=np.zeros((q,q), dtype=complex)
    kList = np.zeros(q)
    cent=np.int(q/2)
    for i in np.arange(q):
        kList[i]=kj-2.0*np.pi*p*(i-cent)/q
        H[i,i]=-2.0*tx*np.cos(kList[i])
        if i==0:
            if q==1:
                H[i,i]+=-2.0*tm*np.cos(km*q)
            else:
                H[i,i+1]=-tm
                H[i,q-1]=-tm*np.exp(1.0j*km*q)
        elif i==q-1:
            H[i,0]=-tm*np.exp(-1.0j*km*q)
            H[i,i-1]=-tm
        else:
            H[i,i+1]=-tm
            H[i,i-1]=-tm

    return H,kList
    
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

qList = np.arange(2,100,1)
fluxList = np.zeros((0))
ElistAll = np.zeros((0))
for qind, q in enumerate(qList):
    pList = np.arange(q-1)+1
    pList = np.append(1,pList[[np.remainder(q,p)!=0 for p in pList]])
    for pind, p in enumerate(pList):
        flux = np.float(p)/np.float(q)
        Egrid,kjList,kmList,mgrid, Vmag = getEigenspectrum(0.1,0.1,p,q, plot=False)
        Elist = Egrid.flatten()
        fluxList = np.append(fluxList,[flux]*Elist.size)
        ElistAll = np.append(ElistAll,Elist)
        
fig = plt.figure()
pan = fig.add_subplot(111)
pan.plot(fluxList,ElistAll,'k.')