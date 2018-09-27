# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 15:07:47 2016

@author: dng5
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 
from scipy import linalg as sLA
from scipy.linalg import block_diag

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in nm
lambdaL = 1064.0e-9 # Lattice wavelenght in nm
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

n=0
S=2
m0=0
k = 0.0
epsilon=0.02
U=4.4
t =220.0e-6*Erecoil/hbar


cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'

    
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
    return Ftot

def RfLatHam(k, omega, delta, epsilon, U, n, S, m0):
    
    Nlat=2*n+1
    Ntot=Nlat*(2*S+1)
    Kinetic=np.zeros(Ntot)

    for i in range(Ntot):

        latI=np.float(np.divide(i,2*S+1)-n)
        Kinetic[i]=(k-2.0*latI)**2.0
    H=np.diag(Kinetic)
    H+=delta*sLA.block_diag(*[Fz(S)]*Nlat)
    H+=((-1.0)**(S+1))*epsilon*sLA.block_diag(*[Fz(S)**2.0]*Nlat)
    H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
        
    for i in range(Ntot):
        for z in range(Ntot):
            if np.abs(i-z)==(2*S+1):
                H[i,z]=U/4.0         
    return H


def rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,rampOnt=0.00000015*Erecoil/hbar,steps=50,adiabatic=False):
    tlist=np.linspace(0.0,rampOnt,steps)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = sLA.eigh(RfLatHam(k, 0.0, 0.0, epsilon, U, n, S, m0))
    if adiabatic:
        EnergyA, VA = sLA.eig(RfLatHam(k, omega,delta, epsilon, U, n, S, m0))
        sort=np.argsort(EnergyA)
        Esorted, eVsorted = EnergyA[sort], VA[:,sort]        
        psi0=eVsorted[:,0]
    else:
#        psi0=np.zeros(Energy1.size).reshape(2*n+1,2*S+1).transpose()
#        psi0[m0+S]=np.ones(2*n+1)/np.sqrt(np.float(2*n+1))
#        psi0=psi0.transpose().flatten()
        
        
        sort=np.argsort(Energy1)
        V1sorted=V1[:,sort]
        psiLat=V1sorted[:,0].reshape(2*n+1,2*S+1).transpose()
        psi0=np.zeros(Energy1.size).reshape(2*n+1,2*S+1).transpose()
        psi0[m0+S]=psiLat[S]/np.sqrt(np.sum(psiLat[S]**2.0))
        psi0=psi0.transpose().flatten()
        if any(psi0 != psi0):
            print 'nan in psi0!'
        
        for t in tlist:
            omegaLoc=omega*t/rampOnt
            Energy,V=sLA.eig(RfLatHam(k, omegaLoc, delta, epsilon, U, n, S, m0))
            V=V+0.0j
            Vinv=np.conj(V.transpose())
            psi0=np.dot(Vinv,psi0)
            teo=np.diag(np.exp(-1.0j*Energy*dt))
            psi0=np.dot(teo,psi0)
            psi0=np.dot(V,psi0)
        
  #  psi0=psi0.reshape(n,3).transpose()
   # pops=np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return psi0

def propagateRfHamiltonian(t, k,omega,delta,epsilon,U,n,S,m0):  

    t=np.array(t)
    if n==0:
        psi0=np.zeros(2*S+1,dtype=complex)
        psi0[m0+S]=1.0+0.0j
    else:
        psi0 = rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0)
    H = RfLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energy, V = sLA.eig(H)

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
    spinPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1),axis=1).flatten() 
    return spinPops

def propagateRfHofDelta(rfFreq,resRfFreq,t,k, omega,epsilon,U,n,S,m0):  
    if S==2:
        delta=(resRfFreq-rfFreq)
    else:
        delta=-(resRfFreq-rfFreq)
    

#    if n==0:
#        psi0=np.zeros(2*S+1,dtype=complex)
#        psi0[m0+S]=1.0+0.0j
#    else:
#        psi0 = rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0)
    pops=np.zeros((delta.size,2*S+1))
    ind=0
    for delt in delta:
#        H = RfLatHam(k, omega, delt ,epsilon,U,n,S,m0)
#        Energy, V = sLA.eig(H)
#        
#        V = V + 1j*0.0
#        Vinv = np.conjugate(np.transpose(V))
#    
#        # np.outer(t, Energy).flatten() creates a matrix for all t
#        Umat = np.diag(np.exp(-1j*t*Energy))  
#    
#        a = np.dot(Vinv, psi0)
#        b = np.dot(Umat, a)                                                            
#        psi = np.dot(V, b)
#        pop=np.absolute(psi)**2.0 
#        spinPop=np.sum(pop.reshape(2*n+1,2*S+1),axis=0)
#        if any(spinPop != spinPop):
#            print 'nan at delta,resRfFreq: '
#            print delt,resRfFreq
#
#        pops[ind]=spinPop
        pops[ind]=propagateRfHamiltonian(t, k,omega,delt,epsilon,U,n,S,m0)
        ind+=1                    
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal

    return pops.flatten()
  
def plotPulsedRFvsDetuning(deltaMax, step, psi0 = [0.0,0.0,1.0,0.0,0.0],U=0.0, k=0.0, omega=4.0, t=0.0, epsilon=0.0):
    deltaList = np.arange(-deltaMax,deltaMax,step)
    pops = np.zeros((deltaList.size,2*S+1))

    t=np.array(t)
    np.append(np.arange(424,448),np.arange(451,454))#
    i=0    
    for delta in deltaList:
        pops[i]=propagateRfHamiltonian(t, k,omega,delta,epsilon,U,n,S,m0)

        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.set_title(r'$\Omega$ = ' + str(omega) + r' $E_r/\hbar$, pulse time = '+str(t)+r' recoils, $\epsilon$ = ' + str(epsilon)+ r' $E_r$')
    for mF in np.arange(-S,S+1):
        pan.plot(deltaList,pops[:,mF+S],cDict[mF], label='mF='+str(mF))

    pan.set_xlabel(r'\delta [$E_r/\hbar$]')
    legend()
    return 
    
  
def plotAdiabaticRFvsDetuning(deltaMax, step,U=0.0, k=0.0, omega=4.0, epsilon=0.0):
    deltaList = np.arange(-deltaMax,deltaMax,step)
    pops = np.zeros((deltaList.size,2*S+1))

    
    i=0    
    for delta in deltaList:
        psi=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,adiabatic=True)
        pop=np.absolute(psi)**2.0                     

        spinPops=np.sum(pop.reshape(2*n+1,2*S+1),axis=0).flatten() 
        pops[i]=spinPops

        i+=1
        
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.set_title(r'$\Omega$ = ' + str(omega) + r' $E_r/\hbar$, $\epsilon$ = ' + str(epsilon)+ r' $E_r$')
    pan.plot(deltaList,pops[:,2*S],'c-', label='mF=+2')
    pan.plot(deltaList,pops[:,2*S-1],'b-', label='mF=+1')
    pan.plot(deltaList,pops[:,2*S-2],'g-', label='mF=0')
    pan.plot(deltaList,pops[:,2*S-3],'r-', label='mF=-1')
    pan.plot(deltaList,pops[:,2*S-4],'m-', label='mF=-2')
    pan.set_xlabel(r'$\delta$ [$E_r/\hbar$]')
    legend()
    return 
    
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
        H=RfLatHam(k, omega, delt, epsilon, U, n, S, m0)
        Energies, eigenstates = sLA.eig(H)
        sort = np.argsort(Energies)
        Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
        E[i]=Esorted
        ev=eVsorted.transpose()
        evA=ev.reshape(N*s,N,s)
        popsA=np.sum(evA*np.conjugate(evA),axis=1)
        pops[i]=popsA
        
#        evB=ev.reshape(N*s,N*s)
#        popsB=evB*np.conjugate(evB)
#        ktot[i]=np.einsum('i,ji->j',kstates,popsB)
        
        m[i]=np.dot(popsA,mFvect)
        i=i+1
   
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title(r'$\Omega$ = '+str(np.round(omega,2))+r'$E_L$, $\epsilon$ = '+str(np.round(epsilon,3))+r'$E_L$, U = '+str(np.round(U,2))+r'$E_L$')
    for i in range(s):   
        d=panel.scatter(deltaList,E[:,i],c=m[:,i],vmin=-S,vmax=S, marker='_')
    panel.set_xlabel(r'$\delta/E_L$')
    panel.set_ylabel(r'$E/E_L$')
    plt.colorbar(d)
#roiList=[[505, 535, 430, 470], [500, 530, 475, 515], [495, 525, 525, 565]]
#
#filestart=4
#filestop=36
#fileroot = 'X:/2016/June/29/PIXIS_29Jun2016' 
#counts,fractions,waveDict = readIgor.batchCountMultipleROIs(fileroot,filestart,filestop,roiList,bgndLoc='top')                      
#rfFreqList=waveDict['rfFreq']
#rfFreqRecoils = rfFreqList*1.0e6*(hbar*2.0*pi)/Erecoil
#fractions=np.array(fractions)
#rfResGuess=0.819*1.0e6*(hbar*2.0*pi)/Erecoil
#a=np.array(rfFreqRecoils)
#
#popt,pcov = optimize.curve_fit(propagateRfHofDelta,a,fractions.flatten(), p0=(rfResGuess,0.9))
#print popt,pcov
#rfForFit=np.linspace(np.min(rfFreqRecoils),np.max(rfFreqRecoils),200)
#pops_fitted=propagateRfHofDelta(rfForFit,*popt)
#pop0 = np.array([pops_fitted[i*3] for i in range(np.int(pops_fitted.size/3))])
#pop1 = np.array([pops_fitted[i*3+1] for i in range(rfForFit.size)]) 
#pop2 = np.array([pops_fitted[i*3+2] for i in range(rfForFit.size)]) 
#
#weightGuess=np.array([26.0/21.0,1.0,26.0/31.0])
#countsW=weightGuess*counts
#total=np.sum(countsW, axis=1)
#fractionsW=np.zeros_like(fractions)
#for i in range(total.size):
#    fractionsW[i] = countsW[i]/total[i]
#
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.set_title('Omega = ' + str(np.round(popt[1],2)) + ' Er/hbar, resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(2.0*pi*hbar),5))+' MHz')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
#panel.plot(rfFreqList,fractions[:,0],'bo', label='mF=+1')
#panel.plot(rfFreqList,fractions[:,1],'go', label='mF=0')
#panel.plot(rfFreqList,fractions[:,2],'ro', label='mF=-1')
#panel.plot(rfForFit*Erecoil*1e-6/(2.0*pi*hbar),pop0,'b-')
#panel.plot(rfForFit*Erecoil*1e-6/(2.0*pi*hbar),pop1,'g-')
#panel.plot(rfForFit*Erecoil*1e-6/(2.0*pi*hbar),pop2,'r-')
#panel.set_xlabel('Rf frequency [MHz]')
datafile=np.load('12Apr2017_F2Calib_newBgnd.npz')        #
rfFreqList=datafile['tlist']
#signalGood=datafile['signalGood']
#imbal=datafile['imbalArray']
#cutoff=0.2
#fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
fieldGoodArray=np.ones(rfFreqList.size,dtype=bool)
#numPoints=datafile['numPoints']
#fieldGoodArray=(numPoints>1)
rfFreqRecoils = rfFreqList[fieldGoodArray]*1.0e6*(hbar*2.0*np.pi)/Erecoil
fractionP=datafile['fractionP'][fieldGoodArray]
fraction0=datafile['fraction0'][fieldGoodArray]
fractionM=datafile['fractionM'][fieldGoodArray]
rfFreqList=rfFreqList[fieldGoodArray]

fractions=np.append(fractionM,np.append(fraction0,fractionP))
#num1=datafile['num1Array'][fieldGoodArray]
#num2=datafile['num2Array'][fieldGoodArray]
#uwaveTot=(num1+num2)/2.0

numP=datafile['numPnorm'][fieldGoodArray]#/uwaveTot
num0=datafile['num0norm'][fieldGoodArray]#/uwaveTot
numM=datafile['numMnorm'][fieldGoodArray]#/uwaveTot
counts=np.append(numM,np.append(num0,numP))
if S==2:
    fractionP2=datafile['fractionP2'][fieldGoodArray]
    fractionM2=datafile['fractionM2'][fieldGoodArray]
    numP2=datafile['numP2norm'][fieldGoodArray]#/uwaveTot
    numM2=datafile['numM2norm'][fieldGoodArray]#/uwaveTot
    fractions=np.append(fractionM2,np.append(fractionM,np.append(fraction0,np.append(fractionP,fractionP2))))
    counts=np.append(numM2,np.append(numM,np.append(num0,np.append(numP,numP2))))
counts=counts.reshape(2*S+1,rfFreqList.size).transpose()
fractions=fractions.reshape(2*S+1,rfFreqList.size).transpose()
rfResGuess=0.5329*1.0e6*(hbar*2.0*np.pi)/Erecoil
a=np.array(rfFreqRecoils)
omega=0.48
def weightRfHofDelta(rfFreq,resRfFreq,wM2,w0,wP,wP2):
    
    weights=np.array([1.0,w0,wP])
    if S==2:
        weights=np.array([wM2,1.0,w0,wP,wP2])
    fracs = np.array(propagateRfHofDelta(rfFreq,resRfFreq,t,k, omega,epsilon,U,n,S,m0))
    fracs = fracs.reshape(rfFreq.size,2*S+1)
    fracsW = (1.0/weights)*fracs
    fracsW=fracsW.transpose()
    fracsW=fracsW/np.sum(fracsW,axis=0)
    fracsW=fracsW.transpose()
    return fracsW.flatten()
    
def weightRfHofDeltaCounts(rfFreq,rfRes,wM2,w0,wP,wP2,numT):#(rfFreq,rfRes,w0,wP,numT)
    #rfRes=rfResGuess
    weights=np.array([1.0,w0,wP])
    if S==2:
        weights=np.array([wM2,1.0,w0,wP,wP2])
    fracs = np.array(propagateRfHofDelta(rfFreq,rfRes,t,k, omega,epsilon,U,n,S,m0))
    count = fracs*numT
    count = count.reshape(rfFreq.size,2*S+1)
    countW = (1.0/weights)*count

    return countW.flatten()

(rfResFit,wM2Fit,w0Fit, wPFit,wP2Fit,numT),pcov = optimize.curve_fit(weightRfHofDeltaCounts,a,counts.flatten(), p0=(rfResGuess,1.0,1.1,1.0,0.9,100.0))
omegaFit=omega
#rfResFit=rfResGuess
print  rfResFit,wM2Fit, w0Fit,wPFit,wP2Fit,numT,np.sqrt(np.diag(pcov))#rfResFit, w0Fit,wPFit,numT,np.sqrt(np.diag(pcov))#
rfForFit=np.linspace(np.min(rfFreqRecoils),np.max(rfFreqRecoils),200)
pops_fitted=propagateRfHofDelta(rfForFit,rfResFit,t,k,omegaFit,epsilon,U,n,S,m0)
pops_fitted=pops_fitted.reshape(rfForFit.size,2*S+1)
#popM = np.array([pops_fitted[i*3] for i in range(np.int(pops_fitted.size/3))])
#pop0 = np.array([pops_fitted[i*3+1] for i in range(rfForFit.size)]) 
#popP = np.array([pops_fitted[i*3+2] for i in range(rfForFit.size)]) 

weightFit=np.array([1,w0Fit,wPFit])
if S==2:
    weightFit=np.array([wM2Fit,1.0,w0Fit,wPFit,wP2Fit])
#weightFit=np.ones(5)
countsW=weightFit*counts
total=np.sum(countsW, axis=1)
fractionsW=np.zeros_like(fractions)
for i in range(total.size):
    fractionsW[i] = countsW[i]/total[i]

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
panel.set_title('Omega = ' + str(np.round(omegaFit,2)) + ' Er/hbar, res at '+str(np.round(rfResFit*Erecoil*1e-6/(2.0*np.pi*hbar),5))+' MHz, w0='+str(np.round(w0Fit,2))+' , wP= '+str(np.round(wPFit,2))+', wM2 = '+str(np.round(wM2Fit,2))+', wP2 = '+str(np.round(wP2Fit,2)))#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
panel.plot(rfFreqList,fractionsW[:,S+1],'bo', label='mF=+1')
panel.plot(rfFreqList,fractionsW[:,S],'go', label='mF=0')
panel.plot(rfFreqList,fractionsW[:,S-1],'ro', label='mF=-1')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pops_fitted[:,S+1],'b-')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pops_fitted[:,S],'g-')
panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pops_fitted[:,S-1],'r-')
panel.plot(rfFreqList,fractions[:,S+1],'bo', fillstyle='none')
panel.plot(rfFreqList,fractions[:,S],'go', fillstyle='none')
panel.plot(rfFreqList,fractions[:,S-1],'ro', fillstyle='none')
if S==2:
    panel.plot(rfFreqList,fractionsW[:,S+2],'co', label='mF=+2')
    panel.plot(rfFreqList,fractionsW[:,S-2],'mo', label='mF=-2')
    panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pops_fitted[:,S+2],'c-')
    panel.plot(rfForFit*Erecoil*1e-6/(2.0*np.pi*hbar),pops_fitted[:,S-2],'m-')
    panel.plot(rfFreqList,fractions[:,S+2],'co', fillstyle='none')
    panel.plot(rfFreqList,fractions[:,S-2],'mo', fillstyle='none')
panel.set_xlabel('Rf frequency [MHz]')
panel.set_ylabel('weighted fractional populations')
legend()



figure=plt.figure()
panel=figure.add_subplot(1,1,1)
#panel.set_title('Omega = ' + str(np.round(popt[1],2)) + ' Er/hbar, resonance at '+str(np.round(popt[0]*Erecoil*1e-6/(2.0*pi*hbar),5))+' MHz')#, epsilon = ' + str(np.round(popt[2],3))+ ' Er')
panel.plot(rfFreqList,counts[:,S+1],'bo',  fillstyle='none')
panel.plot(rfFreqList,counts[:,S],'go',  fillstyle='none')
panel.plot(rfFreqList,counts[:,S-1],'ro',  fillstyle='none')
panel.plot(rfFreqList,np.sum(counts,axis=1),'ko', fillstyle='none')
panel.plot(rfFreqList,countsW[:,S+1],'bo', label='mF=+1')
panel.plot(rfFreqList,countsW[:,S],'go', label='mF=0')
panel.plot(rfFreqList,countsW[:,S-1],'ro', label='mF=-1')
panel.plot(rfFreqList,np.sum(countsW,axis=1),'ko', label='total')
if S==2:
    panel.plot(rfFreqList,counts[:,S+2],'co',  fillstyle='none')
    panel.plot(rfFreqList,counts[:,S-2],'mo',  fillstyle='none')
    panel.plot(rfFreqList,countsW[:,S+2],'co', label='mF=+2')
    panel.plot(rfFreqList,countsW[:,S-2],'mo', label='mF=-2')
panel.set_ylabel('weighted counted atom number')
panel.set_xlabel('Rf frequency [MHz]')
