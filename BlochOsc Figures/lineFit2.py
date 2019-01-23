# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:05:03 2016

@author: dng5
"""
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

def line(x,A,B):
    return A*x+B
    
def lineFit(xvars,yvars,xlabel,ylabel,errorbars=False,yerr=0.0,plot=True,**kwargs):
    if errorbars:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,sigma=yerr,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        if plot:
            figure=plt.figure()
            pan=figure.add_subplot(1,1,1)
            pan.errorbar(xvars,yvars,yerr=yerr,fmt='bo')
            pan.plot(xrangefit,data_fitted,'b-')
            pan.set_title('Fit params in Ax+B, A='+str(np.round(A,3))+'+/-'+str(np.round(dA,4))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
            pan.set_xlabel(xlabel)
            pan.set_ylabel(ylabel)
    else:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        if plot:
            figure=plt.figure()
            pan=figure.add_subplot(1,1,1)
            pan.plot(xvars,yvars,'bo')
            pan.plot(xrangefit,data_fitted,'b-')
            pan.set_title('Fit params in Ax+B, A='+str(np.round(A,3))+'+/-'+str(np.round(dA,4))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
            pan.set_xlabel(xlabel)
            pan.set_ylabel(ylabel)
    return A,B,dA,dB
    
def parabola(x,a,b,c):
    return a*(x**2.0)+b*x+c
    
def lor(x,x0,A,Gamma):
    out=A/((x-x0)**2.0+(Gamma/2.0)**2.0)
    return out
    
def lorFit(xvars,yvars,xlabel,ylabel,p0=(0,1,1)):
    (a,b,c), cov=optimize.curve_fit(lor,xvars,yvars)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=lor(xrangefit,*(a,b,c))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (x0,A,Gamma): '+str(np.round(a,8))+', '+str(np.round(b,5))+', '+str(np.round(c,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return a,b,c

def gauss(x,x0,A,sigma,offset):
    out=A*np.exp(-(x-x0)**2.0/(2.0*sigma**2.0)) + offset
    return out
    
def gaussFit(xvars,yvars,xlabel,ylabel,p0=(0,1,1,0)):
    (x0,A,sigma,offset), cov=optimize.curve_fit(gauss,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=gauss(xrangefit,x0,A,sigma,offset)
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (x0,A,sigma,offset): '+str(np.round(x0,3))+', '+str(np.round(A,3))+', '+str(np.round(sigma,3))+', '+str(np.round(offset,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return x0,A,sigma,offset
    
def parabolicFit(xvars,yvars,xlabel,ylabel,p0=(1,1,0)):
    (a,b,c), cov=optimize.curve_fit(parabola,xvars,yvars)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=parabola(xrangefit,*(a,b,c))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (a,b,c) in a*x^2+b*x+c='+str(np.round(a,8))+', '+str(np.round(b,5))+', '+str(np.round(c,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return a,b,c
    
def sine(x,A,f,phi,offset):
    return offset+A*np.sin(f*x*2.0*np.pi+phi)
    
def sine60(x,A,phi,offset):
    return offset+A*np.sin(60*x*2.0*np.pi+phi)
    
def sineFit(xvars,yvars,xlabel,ylabel,p0=(1,1,0,0),plot=True):
    (A,f,phi,offset), cov=optimize.curve_fit(sine,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=sine(xrangefit,*(A,f,phi,offset))
    yfit=sine(xvars,A,f,phi,offset)
    resids=yvars-yfit
    chi2=np.sum(resids**2.0/yfit)
    uncerts = np.sqrt(np.diag(cov))

    if plot:
        figure=plt.figure()
        pan=figure.add_subplot(1,1,1)
        pan.plot(xvars,yvars,'bo')
        pan.plot(xrangefit,data_fitted,'b-')
        pan.set_title('Fit params (A,f,phi,offset)='+str(np.round(A,3))+', '+str(np.round(f,3))+', '+str(np.round(phi,4))+', '+str(np.round(offset,3)))
        pan.set_xlabel(xlabel)
        pan.set_ylabel(ylabel)
    return A,f,phi,offset,resids,chi2,uncerts
    
def sineFit60(xvars,yvars,xlabel,ylabel,p0=(1,0,0)):
    (A,phi,offset), cov=optimize.curve_fit(sine60,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=sine60(xrangefit,*(A,phi,offset))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (A,phi,offset)='+str(np.round(A,3))+', '+str(np.round(phi,4))+', '+str(np.round(offset,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return A,phi,offset
    
def TwoSine60(x,A,phi,offset):
    return offset+A*np.sin(60*x*2.0*np.pi+phi)+A*np.sin(60*x*2.0*np.pi+phi)
    
def TwoSineFit60(xvars,yvars,xlabel,ylabel,p0=(1,0,0)):
    (A,phi,offset), cov=optimize.curve_fit(sine60,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=sine60(xrangefit,*(A,phi,offset))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (A,phi,offset)='+str(np.round(A,3))+', '+str(np.round(phi,4))+', '+str(np.round(offset,3)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return A,phi,offset
    
def plane((x,y), A,B,C):
    pl=A*x+B*y+C
    return pl.ravel()    
    
def expDecay(t,tau,A,offset):
    return A*np.exp(-t/tau) + offset
    
def expFit(xvars,yvars,xlabel,ylabel,p0=(1,1,0)):
    (tau,A,offset), cov=optimize.curve_fit(expDecay,xvars,yvars,p0=p0)
    xrangefit=np.linspace(np.min(xvars),np.max(xvars),600)
    data_fitted=expDecay(xrangefit,*(tau,A,offset))
    print np.sqrt(np.diag(cov))
    figure=plt.figure()
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.plot(xrangefit,data_fitted,'b-')
    pan.set_title('Fit params (tau,A,offset)='+str(np.round(tau,3))+', '+str(np.round(A,3))+', '+str(np.round(offset,4)))
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return tau,A,offset
    
def easyPlot(xvars,yvars,xlabel='',ylabel='',title=''):
    figure=plt.figure()    
    pan=figure.add_subplot(1,1,1)
    pan.plot(xvars,yvars,'bo')
    pan.set_title(title)
    pan.set_xlabel(xlabel)
    pan.set_ylabel(ylabel)
    return
#
#xCent=np.array([163,161,165,167,165,166,167,162,164,160,159,163,160,159])
#xGuess2=np.array([112,116,123,126,131,135,140,106,101,101,91,87,87,126])
#yCent=np.array([176,182,185,190,194,200,206,170,162,165,154,150,149,198])
##(A,B,C), cov = optimize.curve_fit(plane,(xCent,xGuess2),yCent)
##print A,B,C, np.sqrt(np.diag(cov))
#    
#xrot=np.array([129,139,156,106,138,144])
#y=np.array([116,133,146,94,132,145])
#
#Apd=np.array([7.76,3.4,0.848,0.192])
#Apow=np.array([82.1,36.4,9.1,2.1])
#Cpd=np.array([6.04,3.0,0.776,0.156])
#Cpow=np.array([68.7,34.4,9.0,2.1])
#
#latPD=np.array([1.06,1.26,0.244,0.192,0.292])
#latE=np.array([37.0,43.4,6.82,5.26,8.07])
#
#latPD2=np.array([82,96,104])
#latE2=np.array([4.1,5.2,5.5])
#
#
#ramanApd=np.array([363,363,365,364,365,365,365,364,366,362,362,362,363,362,363,
#                   363,362,362,362,362,362,362,362,358,362,363,362,361,362,362,
#                   361,360,362,362,362,361,362,362,362,362,363,362,362,363,363,
#                   362,362,362,362,362])
#ramanCpd=np.array([377,377,372,373,384,374,378,378,385,376,380,382,381,379,383,
#                   379,384,384,384,380,386,388,385,388,388,388,391,388,389,392,
#                   390,388,384,392,392,392,388,387,385,390,387,390,385,385,387,
#                   385,392,390,387,385])
#
#dataFile=np.load('18Jun2016_files_208-241.npz')
#signalGood=dataFile['signalGood']
#imbal=dataFile['imbalArray']
#
#latCommand=([0.2,0.1,0.06,0.3,0.08])
#latPD=np.array([220.0,118.0,78.0,318.0,98.0])
#latPower=np.array([10.3,5.0,2.5,16.0,3.8])
##lineFit(latPD2,latE2,'Lattice PD [mV]',r'Lattice depth [$E_L$]')
##lineFit(Cpd,Cpow,'RamanC PD [V]',r'RamanC power [mW]')
##cutoff=0.05
##fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
##time=dataFile['tlist'][fieldGoodArray] #time=np.linspace(0,0.005,num=fieldGoodArray.size)[fieldGoodArray]#
##fractionP=dataFile['fractionP'][fieldGoodArray]
##sineFit(time,fractionP,'hold time [s]','fraction in mF=+1',(0.1,60,0.0,0.1))
##fractionM=dataFile['fractionM'][fieldGoodArray]
##sineFit(time,fractionM,'hold time [s]','fraction in mF=-1',(0.1,60,-1.3,0.3))
##fraction0=dataFile['fraction0'][fieldGoodArray]
##sineFit(time,fraction0,'hold time [s]','fraction in mF=0',(0.1,60,-1.3,0.3))
##fig3=plt.figure()
##pan3=fig3.add_subplot(1,1,1)
##pan3.plot(time*1.0e3,fractionP,'bo', label='mF=+1')
##pan3.plot(time*1.0e3,fraction0,'go', label='mF=0')
##pan3.plot(time*1.0e3,fractionM,'ro', label='mF=-1')
##pan3.set_xlabel(r'Oscillation time [ms]')
##pan3.set_ylabel('Spin populations')
###
#filename='16Jun2017_files_304-323.npz'
#dataFile=np.load(filename)
#kick='pos'
#imbal=dataFile['imbalArray']
#signalGood=dataFile['signalGood']
#cutoff=0.5
#fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
#
#time=dataFile['tlist'][fieldGoodArray]
#qlist=dataFile['qlist'][fieldGoodArray]
#sort=np.argsort(time)
#time=time[sort]
#qlist=qlist[sort]
#imbal2=imbal[sort]
#for i in range(qlist.size-1):
#    if kick=='pos':
#        if qlist[i+1]>qlist[i]+0.8:
#            qlist[i+1]=qlist[i+1]-2.0
#    if kick=='neg':
#        if qlist[i+1]<qlist[i]-1.0:
#            qlist[i+1]=qlist[i+1]+2.0
#
#        
##A,B,dA,dB=lineFit(time[:time.size],qlist[:time.size],'odtkick',r'quasimomentum [$k_L$]')
#
#fractionP2=dataFile['fractionP2'][fieldGoodArray][sort]
#fractionP=dataFile['fractionP'][fieldGoodArray][sort]
#fraction0=dataFile['fraction0'][fieldGoodArray][sort]
#fractionM=dataFile['fractionM'][fieldGoodArray][sort]
#fractionM2=dataFile['fractionM2'][fieldGoodArray][sort]
#
#sineFit60(time,fractionP,'HzDelay3','fractionP',p0=(1,0,0))
#sineFit60(time,fractionM,'HzDelay3','fractionM',p0=(1,0,0))
#sineFit60(time,fractionP2,'HzDelay3','fractionP2',p0=(1,0,0))
#sineFit60(time,fractionM2,'HzDelay3','fractionM2',p0=(1,0,0))

#
#ramanCpd=ramanCpd[fieldGoodArray][sort]
#lineFit(ramanCpd,fractionP,'ramanCpd','Fraction in mF=+1')
#lineFit(ramanCpd,fractionM,'ramanCpd','Fraction in mF=-1')
#
#ramanApd=ramanApd[fieldGoodArray][sort]
#lineFit(ramanApd,fractionP,'ramanApd','Fraction in mF=+1')
#lineFit(ramanApd,fractionM,'ramanApd','Fraction in mF=-1')

#print np.std(fractionP2),np.std(fractionP),np.std(fraction0),np.std(fractionM),np.std(fractionM2)
##
#
#figure=plt.figure()
#pan1=figure.add_subplot(1,1,1)
#pan1.plot(ramanCpd,fractionP2,'co', label=r'$m_F$=+2')
#pan1.plot(ramanCpd,fractionP,'bo', label=r'$m_F$=+1')
#pan1.plot(ramanCpd,fraction0,'go', label=r'$m_F$=0')
#pan1.plot(ramanCpd,fractionM,'ro', label=r'$m_F$=-1')
#pan1.plot(ramanCpd,fractionM2,'mo', label=r'$m_F$=-2')
#pan1.set_xlabel('Raman C PD[mV]')
#pan1.set_ylabel('Fractional population')
#legend()
##
#figure=plt.figure()
#pan1=figure.add_subplot(1,1,1)
#pan1.plot(imbal2,fractionP2,'co', label=r'$m_F$=+2')
#pan1.plot(imbal2,fractionP,'bo', label=r'$m_F$=+1')
#pan1.plot(imbal2,fraction0,'go', label=r'$m_F$=0')
#pan1.plot(imbal2,fractionM,'ro', label=r'$m_F$=-1')
#pan1.plot(imbal2,fractionM2,'mo', label=r'$m_F$=-2')
#pan1.set_xlabel('uwave imbalance')
#pan1.set_ylabel('Fractional population')
#legend()
#
#fig2=plt.figure()
#pan2=fig2.add_subplot(1,1,1)
#pan2.plot(qlist-B,2*fractionP2+fractionP-fractionM-2*fractionM2,'bo')
#pan2.set_xlabel('quasimomentum')
#pan2.set_ylabel('magnetization')

#xpos=dataFile['xCenters']
#xpos=xpos[sort]
#fig=plt.figure()
#pan=fig.add_subplot(1,1,1)
#pan.plot(time,xpos,'bo')
#pan.set_title(filename)