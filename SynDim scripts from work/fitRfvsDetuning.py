# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 15:07:47 2016

@author: dng5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 
from numpy import linalg as LA

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in nm
lambdaL = 1064.0e-9 # Lattice wavelenght in nm
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

def RfHamiltonian(k, omega, delta, epsilon):
    H = np.array([[(k)**2.0 - delta, omega/2.0, 0.0],[omega/2.0, k**2.0-epsilon,omega/2.0],[0.0,omega/2.0,(k)**2.0 + delta]])
    return H
    
from scipy.linalg import block_diag



def propagateRfHamiltonian(t, omega, delta):  
    k = 0.0
    epsilon=0.02
    psi0 = np.array([0+1j*0.0, 1.0+1j*0.0, 0.0+1j*0.0])
    H = RfHamiltonian(k, omega, delta ,epsilon)
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
    VV = block_diag(*([V]*t.size))                          
    psi = np.dot(VV, b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    return pops

def propagateRfHofDelta(rfFreq,resRfFreq,omega):  
    k = 0.0
    epsilon=0.02
    t =300.0e-6*Erecoil/hbar
    delta=-(resRfFreq-rfFreq)
    psi0 = np.array([0+1j*0.0, 1.0+1j*0.0, 0.0+1j*0.0])

    pops=np.zeros([delta.size,3])
    i=0
    for delt in delta:
        H = RfHamiltonian(k, omega, delt ,epsilon)
        Energy, V = LA.eig(H)
        
        V = V + 1j*0.0
        Vinv = np.conjugate(np.transpose(V))
    
        # np.outer(t, Energy).flatten() creates a matrix for all t
        U = np.diag(np.exp(-1j*t*Energy))  
    
        a = np.dot(Vinv, psi0)
        b = np.dot(U, a)                                                            
        psi = np.dot(V, b)
    
        pops[i]=np.absolute(psi)**2.0 
        i+=1                    
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal

    return pops.flatten()
    
def plotRfBandStructOfDelta(omega, epsilon, U, n,S, m0,deltaList=np.linspace(-1.0,1.0,600)):
    i=0  
        
    s=2*S+1
    N=2*n+1
    E=np.zeros((deltaList.size,s*N))
    m=np.zeros((deltaList.size,s*N))
    pops=np.zeros((deltaList.size,s*N,s))
    ktot=np.zeros((deltaList.size,s*N))
    
    mFvect=np.array([np.float(j-S) for j in range(s)])
    for delt in deltaList:
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
   
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\delta$ = '+str(np.round(delta,3))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
    for i in range(s):   
        d=panel.scatter(kList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='_')
    panel.set_xlabel(r'$q/k_L$')
    panel.set_ylabel(r'$E/E_L$')
    plt.colorbar(d)
  
def plotPulsedRFvsDetuning(deltaMax, step, psi0 = [0.0,1.0,0.0], k=0.0, omega=4.0, t=0.0, epsilon=0.0):
    deltaList = np.arange(-deltaMax,deltaMax,step)
    pop0 = np.zeros(deltaList.size)
    pop1 = np.zeros(deltaList.size)
    pop2 = np.zeros(deltaList.size)
    t=np.array(t)
    
    i=0    
    for delta in deltaList:
        p0,p1,p2=propagateRfHamiltonian(t,omega,delta)
     
        pop0[i]=p0
        pop1[i]=p1
        pop2[i]=p2
        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.set_title(r'$\Omega$ = ' + str(omega) + r' $E_r/\hbar$, pulse time = '+str(t)+r' recoils, $\epsilon$ = ' + str(epsilon)+ r' $E_r$')
    pan.plot(deltaList,pop0,'b-', label='mF=+1')
    pan.plot(deltaList,pop1,'g-', label='mF=0')
    pan.plot(deltaList,pop2,'r-', label='mF=-1')
    pan.set_xlabel(r'\delta [$E_r/\hbar$]')
    legend()
    return 
#roiList=[[540, 570, 545, 575], [540, 570, 495, 525],[540, 570, 445, 475]]
#
#filestart=64
#filestop=93
#filelist=np.arange(filestart,filestop+1)
#fileroot = 'Y:/Data/2017/March/10/PIXIS_10Mar2017' 
#counts,fractions,waveDict,probes = batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc='top')              
datafile=np.load('10Mar2017_files_65-93.npz')        
rfFreqList=datafile['tlist']
rfFreqRecoils = rfFreqList*1.0e6*(hbar*2.0*pi)/Erecoil
fractionP=datafile['fractionP']
fraction0=datafile['fraction0']
fractionM=datafile['fractionM']
fractions=np.append(fractionM,np.append(fraction0,fractionP))
fractions=fractions.reshape(3,rfFreqList.size).transpose()

filelist=np.arange(65,94)
num1=datafile['num1Array']
num2=datafile['num2Array']
numP=datafile['numP']/(num1+num2)
num0=datafile['num0']/(num1+num2)
numM=datafile['numM']/(num1+num2)
counts=np.append(numM,np.append(num0,numP)).reshape(3,rfFreqList.size).transpose()

rfResGuess=0.5348*1.0e6*(hbar*2.0*np.pi)/Erecoil
a=np.array(rfFreqRecoils)

def weightRfHofDelta(rfFreq,resRfFreq,omega,w0,wP):

   # omega=0.5
    weights=np.array([1.0,w0,wP])
    fracs = np.array(propagateRfHofDelta(rfFreq,resRfFreq,omega))
    fracs = fracs.reshape(rfFreq.size,3)
    fracsW = (1.0/weights)*fracs
    fracsW=fracsW.transpose()
    fracsW=fracsW/np.sum(fracsW,axis=0)
    fracsW=fracsW.transpose()
    return fracsW.flatten()

def weightRfHofDeltaCounts(rfFreq,resRfFreq,w0,wP,numT):


    weights=np.array([1.0,w0,wP])
    fracs = np.array(propagateRfHofDelta(rfFreq,resRfFreq,omega))
    count = fracs*numT
    count = count.reshape(rfFreq.size,3)
    countW = (1.0/weights)*count

    return countW.flatten()
omega=0.5   
(resFreqFit,w0Fit,wPFit,numTFit),pcov = optimize.curve_fit(weightRfHofDeltaCounts,a,counts.flatten(), p0=(rfResGuess,2.0,1.1,0.2))

print resFreqFit,w0Fit,wPFit,numTFit,pcov
omegaFit=omega 
rfForFit=np.linspace(np.min(rfFreqRecoils),np.max(rfFreqRecoils),200)
pops_fitted=propagateRfHofDelta(rfForFit,resFreqFit,omegaFit)
popM = np.array([pops_fitted[i*3] for i in range(np.int(pops_fitted.size/3))])
pop0 = np.array([pops_fitted[i*3+1] for i in range(rfForFit.size)]) 
popP = np.array([pops_fitted[i*3+2] for i in range(rfForFit.size)]) 

weightGuess=np.array([26.0/21.0,1.0,26.0/31.0])
weightFit=np.array([1,w0Fit,wPFit])
countsW=weightFit*counts
total=np.sum(countsW, axis=1)
fractionsW=np.zeros_like(fractions)
for i in range(total.size):
    fractionsW[i] = countsW[i]/total[i]

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
panel.set_title('Omega = ' + str(np.round(omegaFit,2)) + ' Er/hbar, resonance at '+str(np.round(resFreqFit*Erecoil*1e-6/(2.0*np.pi*hbar),5))+' MHz, w0='+str(np.round(w0Fit,2))+' , wP= '+str(np.round(wPFit,2)))#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
panel.plot(rfFreqList,fractionsW[:,2],'bo', label='mF=+1')
panel.plot(rfFreqList,fractionsW[:,1],'go', label='mF=0')
panel.plot(rfFreqList,fractionsW[:,0],'ro', label='mF=-1')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),popP,'b-')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pop0,'g-')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),popM,'r-')
panel.plot(rfFreqList,fractions[:,2],'bo', fillstyle='none')
panel.plot(rfFreqList,fractions[:,1],'go', fillstyle='none')
panel.plot(rfFreqList,fractions[:,0],'ro', fillstyle='none')
panel.set_xlabel('Rf frequency [MHz]')
panel.set_ylabel('weighted fractional populations')
legend(loc=2)

#popsWFit=weightRfHofDelta(rfForFit,resFreqFit,omegaFit, w0Fit,wPFit)
#popsWFit=popsWFit.reshape(rfForFit.size,3)
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
##panel.set_title('Omega = ' + str(np.round(popt[1],2)) + ' Er/hbar, resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(2.0*pi*hbar),5))+' MHz')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
#panel.plot(rfForFit,popsWFit[:,2],'b-', label='mF=+1')
#panel.plot(rfForFit,popsWFit[:,1],'g-', label='mF=0')
#panel.plot(rfForFit,popsWFit[:,0],'r-', label='mF=-1')
#
#panel.set_ylabel('counted atom number')
#panel.set_xlabel('Rf frequency [MHz]')
#legend()

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
#panel.set_title('Omega = ' + str(np.round(popt[1],2)) + ' Er/hbar, resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(2.0*pi*hbar),5))+' MHz')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
panel.plot(rfFreqList,counts[:,2],'bo',  fillstyle='none')
panel.plot(rfFreqList,counts[:,1],'go',  fillstyle='none')
panel.plot(rfFreqList,counts[:,0],'ro',  fillstyle='none')
panel.plot(rfFreqList,np.sum(counts,axis=1),'ko', fillstyle='none')
panel.plot(rfFreqList,countsW[:,2],'bo', label='mF=+1')
panel.plot(rfFreqList,countsW[:,1],'go', label='mF=0')
panel.plot(rfFreqList,countsW[:,0],'ro', label='mF=-1')
panel.plot(rfFreqList,np.sum(countsW,axis=1),'ko', label='total')
panel.set_ylabel('counted atom number/uwave atom number')
panel.set_xlabel('Rf frequency [MHz]')
legend()




figure=plt.figure()
panel=figure.add_subplot(1,1,1)
#panel.set_title('Omega = ' + str(np.round(popt[1],2)) + ' Er/hbar, resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(2.0*pi*hbar),5))+' MHz')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')

panel.plot(filelist,np.sum(counts,axis=1),'ko', label='total')
panel.plot(filelist,num1,'bo', label='1st uwave')
panel.plot(filelist,num2,'go',label='2nd uwave')
panel.plot(filelist,2000*np.sum(counts,axis=1)/(num1+num2),'ro', label='normed total *2000')
panel.set_ylabel('weighted total atom number')
panel.set_xlabel('atom number in mF=0')
panel.legend()

