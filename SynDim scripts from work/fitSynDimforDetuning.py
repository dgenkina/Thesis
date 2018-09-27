# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 17:43:10 2016

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


c=0 #1064.0/790.0#
n=2  #2 or 7#
k=0.0
epsilon=0.02
U=4.4
S=2
m0=0


filename = '23Aug2017_files_94-123'
dataFile=np.load(filename+'.npz')
xAxisLabel='HzDelay3'
omega=0.51
pulseTime=0.0007*Erecoil/hbar
rampOnt=0.0000003*Erecoil/hbar
tau=0.0000003*Erecoil/hbar
adiabatic=False
A1=-0.044
phi1=-0.879


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

    for i in range(Ntot):
        spinI=np.float(i%(2*S+1)-S)
        latI=np.float(np.divide(i,2*S+1)-n)
        Kinetic[i]=(k-2.0*(spinI-m0)*c-2.0*latI)**2.0
    H=np.diag(Kinetic)
    H+=delta*sLA.block_diag(*[Fz(S)]*Nlat)
    H+=((-1.0)**(S+1))*epsilon*sLA.block_diag(*[Fz(S)**2.0]*Nlat)
    H+=(np.sqrt(2.0)/2.0)*omega*sLA.block_diag(*[Fx(S)]*Nlat)    
        
    for i in range(Ntot):
        for j in range(Ntot):
            if np.abs(i-j)==(2*S+1):
                H[i,j]=U/4.0         
    return H


def rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,rampOnt=0.0003*Erecoil/hbar,steps=50):
    tlist=np.linspace(0.0,rampOnt,steps)  
    dt=tlist[1]-tlist[0]
    Energy1, V1 = sLA.eigh(RamanLatHam(k, 0.0, 0.0, epsilon, U, n, S, m0))
    sort=np.argsort(Energy1)
    V1sorted=V1[:,sort]
    if c==0:
        psiLat=V1sorted[:,0].reshape(2*n+1,2*S+1).transpose()
        psi0=np.zeros(Energy1.size).reshape(2*n+1,2*S+1).transpose()
        psi0[m0+S]=psiLat[1]/np.sqrt(np.sum(psiLat[1]**2.0))
        psi0=psi0.transpose().flatten()

    else:
        psi0=V1sorted[:,0]
    
    for t in tlist:
        omegaLoc=expRamp(t,rampOnt,0.0,omega,tau)
        Energy,V=sLA.eigh(RamanLatHam(k, omegaLoc, delta, epsilon, U, n, S, m0))
        V=V+0.0j
        Vinv=np.conj(V.transpose())
        psi0=np.dot(Vinv,psi0)
        teo=np.diag(np.exp(-1.0j*Energy*dt))
        psi0=np.dot(teo,psi0)
        psi0=np.dot(V,psi0)
        
  #  psi0=psi0.reshape(n,3).transpose()
   # pops=np.array([np.dot(psi0[0],np.conj(psi0[0])),np.dot(psi0[1],np.conj(psi0[1])),np.dot(psi0[2],np.conj(psi0[2]))])
    return psi0



    
def propagateRLHamiltonian(t, k, omega, delta, epsilon, U, n,S,m0):  
    t=np.array(t)    
    psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energy, V = sLA.eigh(H)

    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))

    # np.outer(t, Energy).flatten() creates a matrix for all t
    U = np.diag(np.exp(-1j*np.outer(t, Energy).flatten()))  

    a = np.dot(Vinv, psi0)
    # This repeats a so that the shape is consitent with U
    aa = np.outer(np.ones(t.size),a).flatten()
                      
    # Have to add the transpose to make shapesfrom mpl_toolkits.mplot3d import Axes3D match 
    b = np.dot(U, aa)                                     
    # Same block diagonal trick for eigenvector matrix
    VV = sLA.block_diag(*([V]*t.size))                          
    psi = np.dot(VV, b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    latPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1),axis=1).flatten() 
    #populations in each spin state, summed over lattice sites, in time step blocks 
    return spinPops
    
def propagateRLHamiltonian2(t, k, omega, delta, epsilon, U, n,S,m0,**kwargs):  
    t=np.array(t)    
    if n==1:
        psi0=np.zeros(2*S+1)
        psi0[m0+S]=1
    else:
        psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,**kwargs)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energy, V = sLA.eigh(H)
    
    mat=np.array([np.identity(Energy.size)]*t.size)


    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    V=np.array([V]*t.size)
    Vinv=np.array([Vinv]*t.size)
    psi0=np.array([psi0]*t.size)

    aa=np.exp(-1j*np.outer(t, Energy))

    U = np.einsum('ij,ijk->ijk',aa,mat)

    a = np.einsum('ijk,ik->ij',Vinv,psi0)
    b = np.einsum('ijk,ik->ij',U,a)                                                           
    psi = np.einsum('ijk,ik->ij',V,b)
    
    pops=np.absolute(psi)**2.0                     
    # Since you want the first value, need to take every 3rd row 
    # and extract the values you want from the diagonal
    latPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1)[:,np.divide(n,2)-1:np.divide(n,2)+2,:],axis=2).flatten() 
    #populations in the -2k_L, 0, and +2k_L lattice sites, summed over spin sites,in time step blocks
    spinPops=np.sum(pops.reshape(t.size,2*n+1,2*S+1),axis=1).flatten() #[:,0]#
    #populations in each spin state, summed over lattice sites, in time step blocks 
    return spinPops

def propagateRLHamiltonianAllPops(t, k, omega, delta, epsilon, U, n,S,m0,**kwargs):  
    t=np.array(t)    
    if n==1:
        psi0=np.zeros(2*S+1)
        psi0[m0+S]=1
    else:
        psi0=rampedOnPsi(k,omega,delta,epsilon,U,n,S,m0,**kwargs)
#    Energy1, V1 = LA.eig( RamanLatHam(0.0, 0.0, 0.0, epsilon, U, n, S, m0))
#    sort=np.argsort(Energy1)
#    V1sorted=V1[:,sort]
#    psi0=V1sorted[:,0]
#    psi0[np.divide(3*n,2)]=1.0+0.0*1j
    H = RamanLatHam(k, omega, delta, epsilon, U, n, S, m0)
    Energy, V = sLA.eigh(H)
    
    mat=np.array([np.identity(Energy.size)]*t.size)

    V = V + 1j*0.0
    Vinv = np.conjugate(np.transpose(V))
    V=np.array([V]*t.size)
    Vinv=np.array([Vinv]*t.size)
    psi0=np.array([psi0]*t.size)


    U = np.einsum('ij,ijk->ijk',np.exp(-1j*np.outer(t, Energy)),mat)

    a = np.einsum('ijk,ik->ij',Vinv,psi0)
    b = np.einsum('ijk,ik->ij',U,a)                                                           
    psi = np.einsum('ijk,ik->ij',V,b)
    
    pops=np.absolute(psi)**2.0                     

    return pops
    
def plotPulsedPops(tf, k, omega, delta, epsilon, U, n,S=2,m0=0,**kwargs):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    spinpops=propagateRLHamiltonian2(tlist, k, omega, delta, epsilon, U, n,S,m0,**kwargs)
    spinpops=spinpops.reshape(tlist.size,2*S+1).transpose()

    
    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omega,delta,U))
    for s in range(2*S+1):
        panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[s],label=r'$m_F$='+str(int(s-S)))
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[1],'r-',label=r'$m_F$=-1')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[2],'g-',label=r'$m_F$=0')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[3],'b-',label=r'$m_F$=+1')
#    panel.plot(tlist*hbar*1.0e6/Erecoil,spinpops[4],'c-',label=r'$m_F$=+2')
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    return
    
    
def plotPulsedPopsGen(tf, k, omega, delta, epsilon, U, n,orders=[(0,-2),(-1,-2),(1,-2)],S=2,m0=0,**kwargs):
    tlist=np.linspace(0,tf*Erecoil/hbar,num=100)
    popsAll=propagateRLHamiltonianAllPops(tlist, k, omega, delta, epsilon, U, n,S,m0,**kwargs)
    popsAll=popsAll.reshape(tlist.size,2*n+1,2*S+1)
    pops=np.zeros((tlist.size,len(orders)))
    for i in range(len(orders)):
        nn=n-orders[i][0]
        ss=S+orders[i][1]
        pops[:,i]=popsAll[:,nn,ss]

    figure=plt.figure()
    panel=figure.add_subplot(1,1,1)
    panel.set_title( r'$\Omega$ = %.3f $E_L/h$, $\delta$ = %.3f $E_L/h$, U= %.2f $E_L$' %(omega,delta,U))
    for i in range(len(orders)):
        panel.plot(tlist*hbar*1.0e6/Erecoil,pops[:,i],label=r'($k_L$,$m_F$)='+str(orders[i]))
    panel.set_xlabel('Lattice pulse time [us]')
    plt.legend()
    return
    
def adiabaticPops(k,omega,delta,epsilon,U,n,S,m0):
    H=RamanLatHam(k, omega, delta, epsilon, U, n,S,m0)
    Energies, eigenstates = LA.eig(H)
    sort = np.argsort(Energies)
    Esorted, eVsorted = Energies[sort], eigenstates[:,sort]
    ev1=eVsorted[:,0].reshape(2*n+1,2*S+1)
    m1pops=np.einsum('ij,ij->j',ev1,np.conjugate(ev1))#np.array([np.dot(ev1[:,0],ev1[:,0]),np.dot(ev1[:,1],ev1[:,1]),np.dot(ev1[:,2],ev1[:,2])])
    return m1pops
    
def plotAdiabatiPopsOfDelta(deltaMin,deltaMax,deltaStep,k,omega,epsilon,U,n,S,m0,t,rampOnt=0.015*Erecoil/hbar):
    dlist=np.arange(deltaMin,deltaMax,deltaStep)
    fAd=np.zeros((dlist.size,2*S+1))
    fRamped=np.zeros((dlist.size,2*S+1))


    for ind,d in enumerate(dlist):
        pops=adiabaticPops(k,omega,d,epsilon,U,n,S,m0)
        fAd[ind]=pops
        pops2=propagateRLHamiltonian2(t, k, omega, d, epsilon, U, n,S,m0,rampOnt=rampOnt)
        fRamped[ind]=pops2

    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(dlist,fAd[:,i], label='mF='+str(i-S))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Fractional pop, adiabatic')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    for i in range(2*S+1):
        pan.plot(dlist,fRamped[:,i], label='mF='+str(i-S))
    pan.set_xlabel(r'Detuning [$E_r/h$]')
    pan.set_ylabel('Fractional pop, ramped')
    pan.set_title(r'$\Omega$=%.3f, U=%.3f,$\epsilon$=%.2f,k=%.2f,n=%.0f,$\tau_{Raman}$=%.3f,$\phi=$%.2f'%(omega,U,epsilon,k,n,rampOnt*hbar/Erecoil,c))
    legend()
    return
def TwoSine60(x,A,phi,offset):
    return offset+A*np.sin(60*x*2.0*np.pi+phi)+A1*np.sin(60*x*2.0*np.pi+phi1)
    
def TwoSineFit60(xvars,yvars,xlabel,ylabel,p0=(1,0,0)):
    (A,phi,offset), cov=optimize.curve_fit(TwoSine60,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=TwoSine60(xrangefit,*(A,phi,offset))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (A,phi,offset)='+str(np.round(A,3))+', '+str(np.round(phi,4))+', '+str(np.round(offset,3))+' A1,phi1 = ' + str(np.round(A1,3))+' ,'+str(np.round(phi1,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return A,phi,offset

def Harmonics60(x,D,E,F,G,phi,offset):
    sin60=0#A*np.sin(60.0*x*2.0*np.pi+phi)
    sin120=0#B*np.sin(120.0*x*2.0*np.pi+phi)
    sin180=0#C*np.sin(180.0*x*2.0*np.pi+phi)
    sin240=D*np.sin(240.0*x*2.0*np.pi+phi)
    sin300=E*np.sin(300.0*x*2.0*np.pi+phi)
    sin360=F*np.sin(360.0*x*2.0*np.pi+phi)
    sin420=G*np.sin(420.0*x*2.0*np.pi+phi)
    return offset+sin60+sin120+sin180+sin240+sin300+sin360+sin420    
    
def HarmonicsFit(xvars,yvars,xlabel,ylabel,p0=(0.1,0.1,0.1,0.1,0.0,0.0)):
    (D,E,F,G,phi,offset), cov=optimize.curve_fit(Harmonics60,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=Harmonics60(xrangefit,*(D,E,F,G,phi,offset))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (D,E,F,G,phi,offset)=  %.2f, %.2f, %.2f, %.2f, %.3f, %.2f'%(D,E,F,G,phi,offset))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return D,E,F,G,phi,offset
signalGood=dataFile['signalGood']
imbal=dataFile['imbalArray']
cutoff=0.5
fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
#fieldGoodArray=np.ones(dataFile['tlist'].size,dtype=bool)


toffs=0.0
tList=dataFile['tlist'][fieldGoodArray]+toffs
sort=np.argsort(tList)
tList=tList[sort]

fractionP=dataFile['fractionP'][fieldGoodArray][sort]
fraction0=dataFile['fraction0'][fieldGoodArray][sort]
fractionM=dataFile['fractionM'][fieldGoodArray][sort]


maxInd=np.int(fractionP.size)
fractions=np.append(np.append(fractionM[:maxInd],fraction0[:maxInd]),fractionP[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']  
 



if S==2:
    fractionP2=dataFile['fractionP2'][fieldGoodArray][sort]
    fractionM2=dataFile['fractionM2'][fieldGoodArray][sort]
    fractions=np.append(np.append(np.append(np.append(fractionM2[:maxInd],fractionM[:maxInd]),fraction0[:maxInd]),fractionP[:maxInd]),fractionP2[:maxInd])#np.load('SynDimCalibData31Mar16.npz')['fractions']   
 
fractions=fractions.reshape(2*S+1,fraction0[:maxInd].size).transpose()   


#tList=waveDict['pulseDelay']#np.load('SynDimCalibData31Mar16.npz')['tList']
tRecoils = np.array(tList*Erecoil/hbar)
fractions=np.array(fractions)
frac0=np.array(fractions[:,0])
#print tRecoils
#psi0=np.array([0+1j*0.0,1.0+1j*0.0,0.0+1j*0.0])



def HofT(t,delt):
    out = propagateRLHamiltonian2(t, k, omega, delt, epsilon, U,n,S,m0,rampOnt=rampOnt)
    return out

def HofTadiabatic(delt):
    out = adiabaticPops(k,omega,delt,epsilon,U,n,S,m0)
    return out
    
yesList=np.array([0,1,2,3,4]) 
def leastSqFromList(deltaList,fractionsList,fractions):
    fracs=np.outer(np.ones(deltaList.size),fractions)
    resids=fractionsList-fracs
    residsSum=np.sum((resids[:,yesList])**2.0,axis=1)
    ind=np.where(residsSum==np.min(residsSum))
    delta=deltaList[ind]
    residuals=resids[ind]    
    return delta, residuals
    
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'
 
'''Comment out this part to avoid recalculating lookup table'''   
#deltaList=np.arange(-0.4,0.4,0.01)
#fractionsList=np.zeros((deltaList.size,2*S+1))
#for ind,d in enumerate(deltaList):
#    if adiabatic:
#        fractionsList[ind]=HofTadiabatic(d)
#        undetunedFracs=HofTadiabatic(0.0)
#    else:
#        fractionsList[ind]=HofT(pulseTime,d)
#        undetunedFracs=HofT(pulseTime,0.0)
'''End lookup table calculation'''
    
deltaFromData=np.zeros(tList.size)
residualsFromData=np.zeros((tList.size,2*S+1))
fittedFractions=np.zeros((tList.size,2*S+1))
for i in range(tList.size):
    
    deltaFromData[i],residualsFromData[i] = leastSqFromList(deltaList,fractionsList,fractions[i])
    fittedFractions[i]=fractionsList[np.where(deltaList==deltaFromData[i])[0][0]]
    

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
for mF in np.arange(-S,S+1):
    
    panel.plot(tList,residualsFromData[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))

panel.set_xlabel(xAxisLabel)
panel.set_ylabel('Residuals from detuning fit')
legend()

A,phi,offset=sineFit60(tList,deltaFromData,xAxisLabel,'fitted detuning',p0=(1,0,0))
detun60List=sine60(tList,A,phi,offset)

#fittedFractions60=np.zeros((tList.size,2*S+1))
#for i in range(tList.size):
#    fittedFractions60[i]=fractionsList[np.where(np.abs(deltaFromData-detun60List[i])==np.min(np.abs(deltaFromData-detun60List[i])))[0][0]]
A,phi,offset=TwoSineFit60(tList,deltaFromData,xAxisLabel,'fitted detuning',p0=(1,0,0))
#D,E,F,G,phi,offset=HarmonicsFit(tList,deltaFromData,xAxisLabel,'fitted detuning')


#xrangefit=np.linspace(np.min(tList)-0.03,np.max(tList)+0.03,600)
#data_fitted=Harmonics60(xrangefit,*(D,E,F,G,phi,offset))

#figure=plt.figure()
#pan=figure.add_subplot(1,1,1)
#pan.plot(tList,deltaFromData,'bo')
#pan.plot(xrangefit,data_fitted,'b-')


stddev=np.std(deltaFromData)
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.plot(tList,deltaFromData,'bo')
#panel.set_xlabel(xAxisLabel)
#panel.set_ylabel('fitted detuning')

figure=plt.figure()
panel=figure.add_subplot(1,1,1)
if S==2:
    panel.plot(tList,fractionP2,'co')
    panel.plot(tList,fractionM2,'mo')
panel.plot(tList,fractionP,'bo')
panel.plot(tList,fraction0,'go')
panel.plot(tList,fractionM,'ro')
for mF in np.arange(-S,S+1):
    panel.plot(tList,fittedFractions[:,mF+S],cDict[mF])
panel.set_xlabel(xAxisLabel)
panel.set_ylabel('Fractional populations')

#
figure=plt.figure()
panel=figure.add_subplot(1,1,1)
for mF in np.arange(-S,S+1):
    
    panel.plot(deltaList,fractionsList[:,mF+S],cDict[mF], label=r'$m_F$='+str(mF))

panel.set_xlabel('Detuning [E_L]')
panel.set_ylabel('Fractions')
legend()

#np.savez(filename+'_detuningFit', residualsFromData=residualsFromData,deltaFromData=deltaFromData, fractions=fractions,tList=tList,omega=omegaGuess,stddev=stddev, undetunedFracs=undetunedFracs,fittedFractions=fittedFractions,fittedFractions60=fittedFractions60)