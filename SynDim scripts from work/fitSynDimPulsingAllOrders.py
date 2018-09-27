# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 14:25:13 2016

@author: dng5
"""


import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy


def RamanLatHamiltonian(k, omega, delta, epsilon, U, n):
    n=np.int(n)
    if n%2==0:
        print "Number of lattice orders n must be odd!"
        
    c=1.064/0.79
    H=np.zeros((3*n,3*n))
    for i in np.arange(3*n):
        spinI=i%3-1
        latI=np.divide(i,3)-np.divide(n,2)
        for j in np.arange(3*n):
            spinJ=j%3-1
            if i==j:
                spinI=np.float(spinI)
                latI=np.float(latI)
                H[i,j]=(k-2.0*spinI*c-2.0*latI)**2.0+spinI*delta-(1-np.abs(spinI))*epsilon
            if np.abs(i-j)==3:
                H[i,j]=U/4.0
            if ((np.abs(i-j)==1) and (np.abs(spinI-spinJ)==1)):
                H[i,j]=omega/2.0
            
    return H

def rampedOnPsi(k,omega,delta,epsilon,U,n,rampOnt=0.0003*Erecoil/hbar):
    tlist=np.linspace(0.0,rampOnt,50)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = LA.eig(RamanLatHamiltonian(0.0, 0.0, 0.0, 0.0, U, n))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    psi0=V1sorted[:,0]
    
    for t in tlist:
        omegaLoc=omega*t/rampOnt
        Energy,V=LA.eig(RamanLatHamiltonian(k,omegaLoc,delta,epsilon, U, n))
        V=V+0.0j
        Vinv=np.conj(V.transpose())
        psi0=np.dot(Vinv,psi0)
        teo=np.diag(np.exp(-1.0j*Energy*dt))
        psi0=np.dot(teo,psi0)
        psi0=np.dot(V,psi0)
        
  #  psi0=psi0.reshape(n,3).transpose()
   # pops=np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return psi0
    
def propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n):  
    t=np.array(t)    
    psi0=rampedOnPsi(k,omega,delta,epsilon,U,n)
#    Energy1, V1 = LA.eig(RamanLatHamiltonian(0.0, 0.0, 0.0, 0.0, U, n))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
   # psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHamiltonian(k, omega, delta ,epsilon,U,n)
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
    latPops=np.sum(pops.reshape(t.size,n,3)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,n,3),axis=1).flatten() 
    #populations in each spin state, summed over lattice sites, in time step blocks 
    countsMatrixBD=sLA.block_diag(*([countsMatrix])*t.size)
    countsPops=np.array(filter(lambda x: x!=0,np.dot(countsMatrixBD,pops)))
    #populations in each counted roi, as defined by countsList
    return countsPops
    
def plotPulsedPops(tf, k, omega, delta, epsilon, U, n):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=70)
    spinpops=propagateRLHamiltonian(tlist, k, omega, delta, epsilon, U, n)
    pop0 = np.array([spinpops[i*3+1] for i in range(tlist.size)])
    popM = np.array([spinpops[i*3] for i in range(tlist.size)]) 
    popP = np.array([spinpops[i*3+2] for i in range(tlist.size)]) 
    
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/h$, $\delta$ = '+str(np.round(popt[1],3))+r' $E_L/h$, U= '+str(np.round(popt[2],3))+r'$E_L$')
    panel.plot(tlist*hbar*1.0e6/Erecoil,pop0,'g-')
    panel.plot(tlist*hbar*1.0e6/Erecoil,popP,'b-')
    panel.plot(tlist*hbar*1.0e6/Erecoil,popM,'r-')
    panel.set_xlabel('Lattice pulse time [us]')
    return
    

dataFile=np.load('24Aug2016_files_124-183_allOrders.npz')
signalGood=dataFile['signalGood']
imbal=dataFile['imbalArray']
cutoff=0.15
fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
##fieldGoodArray=np.ones(dataFile['tlist'].size,dtype=bool)

tList=dataFile['tlist'][fieldGoodArray]

sort=np.argsort(tList)
tList=tList[sort]
fracAll=dataFile['fracAll'][fieldGoodArray][sort]#[:,:5]

tRecoils = np.array(tList*Erecoil/hbar)

n=7
c=1.064/0.79
k=0.0
epsilon=0.027*c*c
U=6.85
omega=0.6
delta=0
center=np.divide(3*n,2)
countsList=[[np.arange(15)],[center,center+3,center-3,center+6,center-6,center+2,center+5,center-1,center+8,center-4,center-2,center+1,center-5,center+4,center-8]]
countsMatrix=np.zeros((3*n,3*n))
countsMatrix[countsList]=1

def constrainedPropagateHam( k, epsilon, n):
    return lambda t,omega, delta,U: np.array(propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n))
HofT=constrainedPropagateHam( k, epsilon, n)

popt,pcov = optimize.curve_fit(HofT,tRecoils,fracAll.flatten(), p0=(0.61,0.001,5.0))
print popt, np.sqrt(np.diag(pcov))
#
tForFit=np.linspace(np.min(tRecoils),np.max(tRecoils),100)
pops_fitted=HofT(tForFit,*popt)
pops_fitted=pops_fitted.reshape(tForFit.size,len(countsList[1])).transpose()

colorList=('b','g','r','c','m','y','b',(0.5,0.5,0.0),(0.0,0.5,0.5),(0.5,0.0,0.5),(0.5,0.5,0.5),(0.1,0.8,0.1),(0.8,0.1,0.1),(0.1,0.1,0.8),(0.3,0.7,0.1))
figure=plt.figure()
panel=figure.add_subplot(1,1,1)
#panel.set_title( r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/h$, $\delta$ = '+str(np.round(popt[1],3))+r' $E_L/h$, U= '+str(np.round(popt[2],3))+r'$E_L$')
for i in range(len(countsList[1])):
    panel.plot(tForFit*hbar*1.0e3/Erecoil,pops_fitted[i],color=colorList[i])
    panel.plot(tList*1.0e3,fracAll.transpose()[i],'o',color=colorList[i])
panel.set_xlabel('Lattice pulse time [us]')
#legend()