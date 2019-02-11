# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 10:23:49 2016

@author: dng5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as sLA
from scipy import optimize
import time

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy

s=2 #value of F quantum number

a=-40*np.pi
b=40*np.pi
alpha=1.0e-6
xSteps=1024


latticeRampOnt=0.03 # in seconds
ramanRampOnt=0.02# in seconds
omegaMax=0.5
delta0=0.0
deltaN=0.0
phiN=0.2*2.0*np.pi
epsilon=0.02
Vmax=4.4
blochOscTime = 0.0093 #Full 2k_L Bloch oscillation time in seconds
x0max=-hbar/(Erecoil*alpha*blochOscTime) #calculate by: hbar/(Erecoil*alpha*(full Bloch oscillation time)) or ((slope from q linefit)*e3*hbar)/(2.0*Erecoil*alpha)
xHopt=10.0e-6 # in seconds
phi=1064.0/790.0
deltaO=0.0
phiO=phiN+0.1*2.0*np.pi
k0=-0.0

tMax=latticeRampOnt + ramanRampOnt + blochOscTime/2.0# in seconds
tstep=0.1 # in recoils
tau=0.005

dataFile=np.load('..\\Raman\\02Sep2016_FastKick.npz')


def params(t, Vmax, omegaMax, x0Max, delta0):
    latticeRampOntr=latticeRampOnt*Erecoil/hbar
    ramanRampOntr=ramanRampOnt*Erecoil/hbar
    xHoptr=xHopt*Erecoil/hbar
    taur=tau*Erecoil/hbar
    if t<latticeRampOntr:
        omega=0.0
        x0=0.0
        V0=Vmax*t/latticeRampOntr
        updating=True
        delta=delta0+deltaN*np.sin(60*2.0*np.pi*t*hbar/Erecoil+phiN)
    elif t<latticeRampOntr+ramanRampOntr:
        omega=(1.0+deltaO*np.sin(60*2.0*np.pi*t*hbar/Erecoil/60.0+phiO))*omegaMax*(1-np.exp((t-latticeRampOntr)/taur))/(1.0-np.exp(ramanRampOntr/taur))#omegaMax*(t-latticeRampOntr)/ramanRampOntr+
        V0=Vmax
        x0=0
        delta=delta0+deltaN*np.sin(60*2.0*np.pi*t*hbar/Erecoil/60.0+phiN)
        updating=True
    elif t< latticeRampOntr+ramanRampOntr+xHoptr:
        V0=Vmax
        omega=omegaMax*(1.0+deltaO*np.sin(60*2.0*np.pi*t*hbar/Erecoil/60.0+phiO))
        x0=x0Max*(t-latticeRampOntr-ramanRampOntr)/xHoptr
        delta=delta0+deltaN*np.sin(60*2.0*np.pi*t*hbar/Erecoil/60.0+phiN)
        updating=True
    else:
        V0=Vmax
        omega=omegaMax*(1.0+deltaO*np.sin(60*2.0*np.pi*t*hbar/Erecoil/60.0+phiO))
        x0=x0Max
        delta=delta0+deltaN*np.sin(60*2.0*np.pi*t*hbar/Erecoil+phiN)
        if ((deltaN==0.0) & (deltaO==0.0)):
            updating=False
        else:
            updating=True
    return V0,omega,x0,delta,updating   


    
def Vho(x,alpha=1.0,x0=0.0):
    return alpha*((x-x0)**2.0)
    
def vLatHo(x,V0=6.0,alpha=0.00015,x0=0.0):
    return (V0/1.0)*(np.sin(x))**2.0+Vho(x,alpha=alpha,x0=x0)

def Fx(S):
    F=np.zeros((2*S+1,2*S+1))
    for i in range(2*S+1):
        for j in range(2*S+1):
            if np.abs(i-j)==1:
                F[i,j]=(1.0/2.0)*np.sqrt(S*(S+1)-(i-S)*(j-S))
    return F

def Fz(S):
    a=np.arange(np.float(-S),np.float(S+1))
    F=np.diag(a)
    return F
    
def FxUpper(S):
    F=np.zeros((2*S+1,2*S+1))
    for i in range(2*S+1):
        for j in range(2*S+1):
            if j-i==1:
                F[i,j]=(1.0/2.0)*np.sqrt(S*(S+1)-(i-S)*(j-S))
    return F
    
        
def vRaman(x,omega=1.0,delta=0.0,epsilon=0.048,phi=4.0/3.0):
    x=np.array(x)    

    v=np.einsum('i,jk->ijk',omega*np.exp(1.0j*2.0*phi*x),FxUpper(s)*np.sqrt(2.0)/2.0)
    v=v+np.conjugate(v).swapaxes(1,2)#np.array([np.triu(v[i])+np.conjugate(np.triu(v[i],1)).transpose() for i in range(x.size)])#np.array(map(np.triu,v)+np.conjugate(map(np.triu,v,np.ones(x.size))))#
    v+=(delta*Fz(s)+((-1.0)**(s+1))*epsilon*(Fz(s)**2.0))*mat
    return v


def getEigenHam2(a,b,xSteps,V,omega=1.0,delta=0.0,epsilon=0.048,phi=4.0/3.0,*args,**kwargs):
    s=1
    xgrid=np.linspace(a,b,xSteps)
  #  print '1. Xgrid:' +str(xgrid.shape)
    eigE=np.zeros((xgrid.size,2*s+1),dtype=complex)
    eigV=np.zeros((xgrid.size,2*s+1,2*s+1),dtype=complex)
    eigVdagger=np.zeros((xgrid.size,2*s+1,2*s+1),dtype=complex)
    for i,x in enumerate(xgrid):
        Vgrid=np.diag(np.array([V(x,*args,**kwargs)]*(2*s+1)))
        Vspin=vRaman(x,omega=omega,delta=delta,epsilon=epsilon,phi=phi)
        Vtot=Vgrid+Vspin
        eigE[i],eigV[i]=sLA.eigh(Vtot)
        eigVdagger[i]=np.conj(eigV[i]).transpose()

    return eigE,eigV,eigVdagger 


def getEigenHam3(a,b,xSteps,V,omega=1.0,delta=0.0,epsilon=0.048,phi=4.0/3.0,*args,**kwargs):
    s=1
 #   xgrid=np.linspace(a,b,xSteps)

  #  print '1. Xgrid:' +str(xgrid.shape)
    eigE=np.zeros((xgrid.size,2*s+1),dtype=complex)
    eigV=np.zeros((xgrid.size,2*s+1,2*s+1),dtype=complex)
    eigVdagger=np.zeros((xgrid.size,2*s+1,2*s+1),dtype=complex)
    
    Vgrid=np.einsum('i,ijk->ijk',V(xgrid,*args,**kwargs),mat)
    Vspin=vRaman(xgrid,omega=omega,delta=delta,epsilon=epsilon,phi=phi)
    Vtot=Vgrid+Vspin
    eigE,eigV=np.linalg.eig(Vtot)
    eigVdagger=np.swapaxes(np.conj(eigV),1,2)

    return eigE,eigV,eigVdagger 

def splitStepPropagatorEigB2(psi,dt,a,b,eigE,eigV,eigVdagger):
 #   xgrid=np.linspace(a,b,psi.shape[0])
    U=np.exp(-1.0j*eigE*dt/2.0)
    psi1=U*psi
    psi1=np.einsum('ijk,ik->ij',eigV,psi1)
    kgrid=np.fft.fftfreq(xgrid.size,d=((b-a)/xgrid.size/(2.0*np.pi)))
    fft1=np.fft.fft(psi1.transpose())
    psi2=np.exp(-1.0j*dt*(kgrid**2.0))*fft1
    fft2=np.fft.ifft(psi2).transpose()
    psi3=np.einsum('ijk,ik->ij',eigVdagger,fft2)
    psi3=U*psi3
    return psi3
    
def splitStepPropagatorUncoupledSpins(psi,V,dt,a,b,*args,**kwargs):
 #   xgrid=np.linspace(a,b,psi.shape[0])
    Vgrid=V(xgrid,*args,**kwargs)
    U=np.array([np.exp(-1.0j*Vgrid*dt/2.0)]*(2*s+1)).transpose()
    psi1=U*psi
    kgrid=np.fft.fftfreq(xgrid.size,d=((b-a)/xgrid.size/(2.0*np.pi)))
    psi2=np.exp(-1.0j*dt*(kgrid**2.0))*np.fft.fft(psi1.transpose())
    psi3=U*np.fft.ifft(psi2).transpose()
    return psi3

def propagateInTime(psi0,V,a,b,tf,dt,omegaMax=1.0,delta0=0.0,epsilon=0.048,phi=4.0/3.0,x0max=0.0,Vmax=6.0,**kwargs):
  #  xgrid=np.linspace(a,b,psi0.shape[0])
 #   dx=(b-a)/(psi0.shape[0])
    tgrid=np.arange(dt,tf,dt)
    fracs=np.zeros((tgrid.size,2*s+1))
    com0=np.zeros(tgrid.size)
    
    for ind,t in enumerate(tgrid):
        psiMag=psi0*np.conj(psi0)
        fracs[ind] = dx*np.sum(psiMag,axis=0)
        com0[ind]=np.sum(xgrid*(np.sum(psiMag,axis=1)*dx))
 
        V0,omega,x0,delta,updating=params(t,Vmax,omegaMax,x0max,delta0)
        change=0
        if omega==0.0:
            psi0=splitStepPropagatorUncoupledSpins(psi0,V,dt,a,b,V0=V0,x0=x0,**kwargs)
        elif updating:    
            eigE,eigV,eigVdagger=getEigenHam3(a,b,psi0.shape[0],V,omega=omega,delta=delta,epsilon=epsilon,phi=phi,x0=x0,V0=V0,**kwargs)
            psi0eigB=np.einsum('ijk,ik->ij',eigVdagger,psi0)
            psi0eigB=splitStepPropagatorEigB2(psi0eigB,dt,a,b,eigE,eigV,eigVdagger)
            psi0=np.einsum('ijk,ik->ij',eigV,psi0eigB)
        else:
            if change==0:
                eigE,eigV,eigVdagger=getEigenHam3(a,b,psi0.shape[0],V,omega=omega,delta=delta,epsilon=epsilon,phi=phi,x0=x0,V0=V0,**kwargs)
                change+=1
            psi0eigB=np.einsum('ijk,ik->ij',eigVdagger,psi0)
            psi0eigB=splitStepPropagatorEigB2(psi0eigB,dt,a,b,eigE,eigV,eigVdagger)
            psi0=np.einsum('ijk,ik->ij',eigV,psi0eigB)
            
        
    return fracs,com0,tgrid,psi0   
    


xgrid=np.linspace(a,b,xSteps)
kgrid=np.fft.fftfreq(xgrid.size,d=((b-a)/xgrid.size/(2.0*np.pi)))
dx=xgrid[1]-xgrid[0]
psiAn=np.exp(-np.sqrt(alpha)*xgrid**2.0/2.0)
psiAn=psiAn/np.sqrt(np.dot(psiAn,np.conj(psiAn)*dx))
Vgrid=vLatHo(xgrid,alpha=alpha)
mat=np.array([np.identity(2*s+1)]*xgrid.size)

prekick=np.exp(1.0j*k0*xgrid)
psi0=np.swapaxes(np.array([np.zeros(psiAn.size),prekick*psiAn,np.zeros(psiAn.size)],dtype=complex),0,1)
if s==2:
    psi0=np.swapaxes(np.array([np.zeros(psiAn.size),np.zeros(psiAn.size),prekick*psiAn,np.zeros(psiAn.size),np.zeros(psiAn.size)],dtype=complex),0,1)
#fig2=plt.figure()
#pan2=fig2.add_subplot(1,1,1)
#pan2.plot(xgrid,Vgrid)

Vgrid=vLatHo(xgrid, alpha=alpha,x0=19640,V0=0.0)
print np.gradient(Vgrid)[xSteps/2]/(xgrid[1]-xgrid[0])


t1=time.clock()
fracs,com0,tgrid,psiOut=propagateInTime(psi0,vLatHo,a,b,tMax*Erecoil/hbar,tstep,omegaMax=omegaMax,delta0=delta0,epsilon=epsilon,phi=phi,alpha=alpha,Vmax=Vmax,x0max=x0max)
t2=time.clock()
print 'Time propagation completed in %f seconds' %(t2-t1)
psiOut=psiOut.reshape(xgrid.size,2*s+1).transpose()

if s==1:
    fracM=fracs[:,2]
    frac0=fracs[:,1]
    fracP=fracs[:,0]
if s==2:
    fracM2=fracs[:,4]
    fracM=fracs[:,3]
    frac0=fracs[:,2]
    fracP=fracs[:,1]
    fracP2=fracs[:,0]
    
#    
#deltaList = np.arange(-0.15,0.15,0.003)
#fracMofD=np.zeros(deltaList.size)
#frac0ofD=np.zeros(deltaList.size)
#fracPofD=np.zeros(deltaList.size)
#
#
#for ind,delta in enumerate(deltaList):
#    
#    t1=time.clock()
#    fracM,frac0,fracP,com0,tgrid,psiOut=propagateInTime(psi0,vLatHo,a,b,tMax*Erecoil/hbar,tstep,omegaMax=omegaMax,delta0=delta0,epsilon=epsilon,phi=4.0/3.0,alpha=alpha,Vmax=Vmax,x0max=x0max)
#    fracMofD[ind]=fracM[-1]
#    frac0ofD[ind]=frac0[-1]
#    fracPofD[ind]=fracP[-1]
#    t2=time.clock()
#    print 'Time propagation completed in %f seconds for index %i' %(t2-t1,ind)
#    psiOut=psiOut.reshape(xgrid.size,2*s+1).transpose()
#
#
#fig1=plt.figure()
#pan1=fig1.add_subplot(1,1,1)
#pan1.plot(deltaList,fracPofD,'b-', label='mF=+1')   
#pan1.plot(deltaList,frac0ofD,'g-', label='mF=0') 
#pan1.plot(deltaList,fracMofD,'r-', label='mF=-1')
#pan1.set_xlabel('Detuning')
#pan1.set_title(r'$\Omega$=%.2f,V=%.2f,$\delta$=%.3f,$\alpha$=%.0f e-6,$x_0$=%.0f,'%(omegaMax,Vmax,delta,alpha*1e6,x0max)+'\n'+ r'$t_{latramp}$=%.3f,$t_{ramanramp}$=%.3f,xSteps=%.0f,tstep=' %(latticeRampOnt,ramanRampOnt,xSteps)+str(np.round(tstep*1e6*hbar/Erecoil,3)))

#fig=plt.figure()
#pan=fig.add_subplot(1,1,1)
#pan.plot(xgrid,psiOut[0]*np.conj(psiOut[0]),'b-', label='mF=+1')   
#pan.plot(xgrid,psiOut[1]*np.conj(psiOut[1]),'g-', label='mF=0') 
#pan.plot(xgrid,psiOut[2]*np.conj(psiOut[2]),'r-', label='mF=-1')
#pan.plot(xgrid, psiAn*np.conj(psiAn),'k-')
#fig.show()
#
#psiFFt=np.fft.fft(psiOut)
#fig=plt.figure()
#pan=fig.add_subplot(1,1,1)
#pan.plot(kgrid,psiFFt[0]*np.conj(psiFFt[0]),'bo', label='mF=+1')   
#pan.plot(kgrid,psiFFt[1]*np.conj(psiFFt[1]),'go', label='mF=0') 
#pan.plot(kgrid,psiFFt[2]*np.conj(psiFFt[2]),'ro', label='mF=-1')
##


#imbal=dataFile['imbalArray']
#signalGood=dataFile['signalGood']
#cutoff=0.35
#fieldGoodArray=((imbal<cutoff) & signalGood)
#fractionP=dataFile['fractionP']#[fieldGoodArray]
#fraction0=dataFile['fraction0']#[fieldGoodArray]
#fractionM=dataFile['fractionM']#[fieldGoodArray]
#time=dataFile['tlist']+(latticeRampOnt+ramanRampOnt)*hbar/Erecoil

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.plot(tgrid*hbar*1e3/Erecoil,fracP,'b-', label='mF=+1')   
pan1.plot(tgrid*hbar*1e3/Erecoil,frac0,'g-', label='mF=0') 
pan1.plot(tgrid*hbar*1e3/Erecoil,fracM,'r-', label='mF=-1')
if s==2:
   pan1.plot(tgrid*hbar*1e3/Erecoil,fracP2,'c-', label='mF=+1')
   pan1.plot(tgrid*hbar*1e3/Erecoil,fracM2,'m-', label='mF=-1') 
pan1.set_title(r'$\Omega$=%.2f,V=%.2f,$\delta_0$=%.3f,$\delta_N$=%.3f,$\phi_N$=%.2f,$\Delta_{\Omega}$=%.1f %%,$\phi_{rel}$=%.2f,'%(omegaMax,Vmax,delta0,deltaN,phiN,deltaO*100.0,phiO-phi)+'\n'+ r'$t_{lat}$=%.3f,$t_{raman}$=%.3f,$x_{Steps}$=%.0f,$t_{step}$=%.1f us,$\alpha$=%.0f e-6,$x_0$=%.0f,' %(latticeRampOnt,ramanRampOnt,xSteps,tstep*1e6*hbar/Erecoil,alpha*1e6,x0max))
pan1.set_xlim(latticeRampOnt*1e3+ramanRampOnt*1e3,tMax*1e3)
#pan1.plot(time*1.0e3+latticeRampOnt*1e3+ramanRampOnt*1e3,fractionP,'bo', label=r'$m_F$=+1')
#pan1.plot(time*1.0e3+latticeRampOnt*1e3+ramanRampOnt*1e3,fraction0,'go', label=r'$m_F$=0')
#pan1.plot(time*1.0e3+latticeRampOnt*1e3+ramanRampOnt*1e3,fractionM,'ro', label=r'$m_F$=-1')

fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.plot(tgrid*hbar*1e3/Erecoil,com0,'b-')

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
mag = fracP-fracM
if s==2:
    mag = mag + 2.0*fracP2-2.0*fracM2
pan1.plot(tgrid*hbar*1e3/Erecoil,mag,'b-', label='mF=+1')   
pan1.set_title(r'$\Omega$=%.2f,V=%.2f,$\delta_0$=%.3f,$\delta_N$=%.3f,$\phi_N$=%.2f,$\Delta_{\Omega}$=%.1f %%,$\phi_{rel}$=%.2f,'%(omegaMax,Vmax,delta0,deltaN,phiN,deltaO*100.0,phiO-phi)+'\n'+ r'$t_{lat}$=%.3f,$t_{raman}$=%.3f,$x_{Steps}$=%.0f,$t_{step}$=%.1f us,$\alpha$=%.0f e-6,$x_0$=%.0f,' %(latticeRampOnt,ramanRampOnt,xSteps,tstep*1e6*hbar/Erecoil,alpha*1e6,x0max))
pan1.set_xlim(latticeRampOnt*1e3+ramanRampOnt*1e3,tMax*1e3)
pan1.set_ylabel('Magnetization')

def gaussian(xlist,sigma,x0):
    return np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))/(sigma*np.sqrt(2.0*np.pi))

def getMode(S,fracList,sigmaList,centGuess='nan',plot=False):
    pad=1
    mFlist=np.arange(2*S+1)-S
    mF=np.lib.pad(mFlist, (pad,pad), 'reflect', reflect_type='odd')
    frac=np.lib.pad(fracList, (pad,pad),'constant',constant_values=(0,0))
    yerr=np.lib.pad(sigmaList, (pad,pad),'constant',constant_values=(0.02,0.02))
    
    if centGuess=='nan':
        centGuess=float(mF[np.where(np.max(frac)==frac)[0][0]])
    
    (sigma,x0), cov=optimize.curve_fit(gaussian,mF,frac,p0=(S/(S+1.0),centGuess),sigma=yerr,absolute_sigma=False)
    (dSigma,dx0)=np.sqrt(np.diag(cov))
    
    if plot:
        print centGuess
        mFforplot=np.linspace(mF[0],mF[-1],num=100)
        fig=plt.figure()
        pan=fig.add_subplot(111)
        pan.errorbar(mF,frac,yerr=yerr,fmt='bo')
        pan.plot(mFforplot,gaussian(mFforplot,sigma,x0))
    return x0, dx0

modeList = np.zeros(fracP.size)
modeSigmaList = np.zeros(fracP.size)
for ind in range(fracP.size):
    modeList[ind], modeSigmaList[ind] = getMode(s,fracs[ind][::-1],0.02*np.ones(2*s+1),centGuess='nan',plot=False)

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.plot(tgrid*hbar*1e3/Erecoil,modeList,'b-')   
pan1.set_title(r'$\Omega$=%.2f,V=%.2f,$\delta_0$=%.3f,$\delta_N$=%.3f,$\phi_N$=%.2f,$\Delta_{\Omega}$=%.1f %%,$\phi_{rel}$=%.2f,'%(omegaMax,Vmax,delta0,deltaN,phiN,deltaO*100.0,phiO-phi)+'\n'+ r'$t_{lat}$=%.3f,$t_{raman}$=%.3f,$x_{Steps}$=%.0f,$t_{step}$=%.1f us,$\alpha$=%.0f e-6,$x_0$=%.0f,' %(latticeRampOnt,ramanRampOnt,xSteps,tstep*1e6*hbar/Erecoil,alpha*1e6,x0max))
pan1.set_xlim(latticeRampOnt*1e3+ramanRampOnt*1e3,tMax*1e3)
pan1.set_ylabel('Modal position')

np.savez('TDSEkickF2experparams',omegaMax=omegaMax,Vmax=Vmax,latticeRampOnt=latticeRampOnt,
         ramanRampOnt=ramanRampOnt,delta0=delta0,epsilon=epsilon,blochOscTime=blochOscTime,
         x0max=x0max,phi=phi, modeList=modeList,modeSigmaList=modeSigmaList,
         fracs=fracs,com0=com0,tgrid=tgrid,psiOut=psiOut)

