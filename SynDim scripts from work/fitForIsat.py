# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 11:08:06 2016

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize



roi=np.array([580,610,155,185])#np.array([510,550,230,260])#np.array([475,515,265,295])#np.array([545,575,195,225])#np.array([455,470,310,320])#
roiB=np.array([550,590,230,260])
#
roiFlea=([100,135,210,245])#([205,260,365,400])

fileroot = 'Y:/Data/2017/March/24/PIXIS_24Mar2017'  
filerootFlea='Y:/Data/2017/March/24/Flea3_24Mar2017'  
filelist=np.arange(64,94)#np.append(np.arange(34,43),np.arange(44,64))#
key='PicIntensity'

def zeroNansInfs(x):
    if x!=x:
        return 0
    elif np.isinf(x):
        return 0
    else:
        return x
zeroNansInfsVector = np.vectorize(zeroNansInfs, otypes=[np.float])

def getRoi(array,xCent,yCent,r=16,eps=5):
    ylen=array.shape[0]
    xlen=array.shape[1]
#    bbox=(xCent-r,yCent-np.int(np.sqrt(r*r-eps*eps)),xCent+r,yCent+np.int(np.sqrt(r*r-eps*eps)))
    
    y,x = np.ogrid[-yCent:ylen-yCent, -xCent:xlen-xCent]
    mask = np.sqrt((x-eps)**2.0 + y*y) +np.sqrt((x+eps)**2.0 + y*y) <= 2.0*r
    maskL = ((mask) & (x<0))
    maskR = ((mask) & (x>0))
    maskT=((mask) & (y<0))
    maskB=((mask) & (y>0))

    counts=np.sum(array[mask])
    countsL=np.sum(array[maskL])
    countsR=np.sum(array[maskR])
    countsT=np.sum(array[maskT])
    countsB=np.sum(array[maskB])

    return (counts, countsL, countsR, countsT, countsB)
    

def findBox(od,xGuess,yGuess):
    allcounts=getRoi(od, xGuess,yGuess)
#    print allcounts
    i=0
    while ((np.abs(allcounts[4]-allcounts[3])>1.1) & (i<20)):
        if (allcounts[4]-allcounts[3])>0:
            yGuess=yGuess+1
#            print "new yGuess = " +str(yGuess)
            allcounts=getRoi(od, xGuess,yGuess)
        else: 
            yGuess=yGuess-1
#            print "new yGuess = " +str(yGuess)
            allcounts=getRoi(od, xGuess,yGuess)
        i=i+1
#    print i
    i=0
    while ((np.abs(allcounts[2]-allcounts[1])>1.1) & (i<20)):
        if (allcounts[2]-allcounts[1])>0:
            xGuess=xGuess+1
#           print "new xGuess = " +str(xGuess)
            allcounts=getRoi(od, xGuess,yGuess)
        else: 
            xGuess=xGuess-1
#           print "new xGuess = " +str(xGuess)
            allcounts=getRoi(od, xGuess,yGuess)
        i=i+1
#    print i
    box=np.array([yGuess-1,yGuess+2,xGuess-1,xGuess+2])
#    print box
    return box
    
def makeImage(fileroot,filenum,roi):
    filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=-41)
    print filename
    od=dict1['OptDepth']
    fig=plt.figure()
    pan=fig.add_subplot(1,3,1)
    pan.imshow(od,vmin=-0.15,vmax=0.5)
    odRoi=od[roi[0]:roi[1],roi[2]:roi[3]]
    pan2=fig.add_subplot(1,3,2)
    pan2.imshow(odRoi,vmin=-0.15,vmax=0.5)
    
    xGuess=(roi[3]-roi[2])/2
    yGuess=(roi[1]-roi[0])/2
    box=findBox(odRoi,xGuess,yGuess)
    pan3=fig.add_subplot(1,3,3)
    pan3.imshow(odRoi[box[0]:box[1],box[2]:box[3]],vmin=-0.15,vmax=0.5)
    return



    
def getIntensities(fileroot,filelist,roi,key,getBox=True):
    PicIntensity=np.zeros(filelist.size)
    Raw1Avg=np.zeros(filelist.size)
    Raw2Avg=np.zeros(filelist.size)
    odAvg=np.zeros(filelist.size)
    for ind,filenum in enumerate(filelist):
        filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
        dict1 =processIBW(filename, angle=-41)
        print filename
        Raw1=dict1['Raw1'][roi[0]:roi[1],roi[2]:roi[3]]
        Raw2=dict1['Raw2'][roi[0]:roi[1],roi[2]:roi[3]]
        od=dict1['OptDepth']
        odRoi=od[roi[0]:roi[1],roi[2]:roi[3]]
        xGuess=(roi[3]-roi[2])/2
        yGuess=(roi[1]-roi[0])/2
        
        if getBox:
            box=findBox(odRoi,xGuess,yGuess)
        else:
            box=roi
        
        
        Raw1Avg[ind]=np.average(zeroNansInfsVector(Raw1[box[0]:box[1],box[2]:box[3]]))
        if Raw1Avg[ind]!=Raw1Avg[ind]:
            Raw1Avg[ind]=0
            print 'Raw1Avg is nan in file number ' + str(filenum)
        Raw2Avg[ind]=np.average(zeroNansInfsVector(Raw2[box[0]:box[1],box[2]:box[3]]))
        if Raw2Avg[ind]!=Raw2Avg[ind]:
            Raw2Avg[ind]=0
            print 'Raw2Avg is nan in file number ' + str(filenum)
        odAvg[ind]=np.average(zeroNansInfsVector(odRoi[box[0]:box[1],box[2]:box[3]]))
        if odAvg[ind]!=odAvg[ind]:
            odAvg[ind]=0
            print 'odAvg is nan in file number ' + str(filenum)
        
        infoString=dict1["Note"]
        waveDictLocal=getIndexedWaves(infoString)
        PicIntensity[ind]=waveDictLocal[key]
        
    return Raw1Avg,Raw2Avg,odAvg,PicIntensity

def isatOptThin(I0, nSigma0, Isat):
    I0minusIf=nSigma0*I0/(Isat+I0)
    return I0minusIf
def line(x,A,B):
    return A*x+B
    
def lineFit(xvars,yvars,xlabel,ylabel,errorbars=False,yerr=0.0,**kwargs):
    if errorbars:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,sigma=yerr,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        figure=plt.figure()
        pan=figure.add_subplot(1,1,1)
        pan.errorbar(xvars,yvars,yerr=yerr,fmt='bo')
        pan.plot(xrangefit,data_fitted,'b-')
        pan.set_title('Fit params in Ax+B, A='+str(np.round(A,12))+'+/-'+str(np.round(dA,12))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
        pan.set_xlabel(xlabel)
        pan.set_ylabel(ylabel)
    else:
        (A,B), cov=optimize.curve_fit(line,xvars,yvars,**kwargs)
        (dA,dB)=np.sqrt(np.diag(cov))
        xrangefit=np.linspace(np.min(xvars),np.max(xvars))
        data_fitted=line(xrangefit,*(A,B))
        
        figure=plt.figure()
        pan=figure.add_subplot(1,1,1)
        pan.plot(xvars,yvars,'bo')
        pan.plot(xrangefit,data_fitted,'b-')
        pan.set_title('Fit params in Ax+B, A='+str(np.round(A,12))+'+/-'+str(np.round(dA,12))+', B='+str(np.round(B,3))+'+/-'+str(np.round(dB,4)))
        pan.set_xlabel(xlabel)
        pan.set_ylabel(ylabel)
    return A,B,dA,dB

date=fileroot.split("/")[-1].split("_")[1]
saveName=date+'_files_'+np.str(filelist[0])+'-'+str(filelist[-1])   
Raw1AvgB,Raw2AvgB,odAvgB,PicIntensity=getIntensities(fileroot,filelist,roiB,key)
Raw1Avg,Raw2Avg,odAvg,PicIntensity=getIntensities(fileroot,filelist,roi,key)
Raw2Avg=Raw2Avg*Raw1AvgB/Raw2AvgB
for i in range(Raw2Avg.size):
    if Raw2Avg[i]==0:
        odAvg[i]=0
    else:
        odAvg[i]=-np.log(Raw1Avg[i]/Raw2Avg[i])

sort=np.argsort(PicIntensity)
minInd=7


(A,B,dA,dB) = lineFit((Raw2Avg[sort][minInd:]-Raw1Avg[sort][minInd:]),odAvg[sort][minInd:],r'$I_0$-$I_f$ [counts]',r'od=-ln($I_f$/$I_0$)')


fig=plt.figure()
fig.suptitle( r'$I_{sat}\alpha^*$ = %.2f+/-%.2f, $\sigma_0 n/\alpha^*$ = %.3f+/-%.3f'%(-1.0/A,dA/(A**2.0),B,dB),size=20 )
pan=fig.add_subplot(1,2,1)
pan.plot(PicIntensity,Raw1Avg,'bo',label=r'$I_f$')
pan.plot(PicIntensity,Raw2Avg,'go',label=r'$I_0$')
pan.plot(PicIntensity,Raw2Avg-Raw1Avg,'ro',label=r'$I_0-I_f$')
pan.set_xlabel('PicIntensity [V]')
pan.set_ylabel('Intensity [counts]')
legend(loc=2)
pan2=fig.add_subplot(1,2,2)
pan2.plot(PicIntensity,odAvg,'bo',label=r'od')
pan2.set_xlabel('PicIntensity [V]')
pan2.set_ylabel('od')


(nSigma0Isat,Isat),cov=optimize.curve_fit(isatOptThin,Raw2Avg[sort][minInd:],Raw2Avg[sort][minInd:]-Raw1Avg[sort][minInd:],p0=(8000,25000.0))
(dnSigma0Isat,dIsat)=np.sqrt(np.diag(cov))
dnSigma0=np.sqrt(((Isat**2.0)*(dnSigma0Isat**2.0)+(nSigma0Isat**2.0)*(dIsat**2.0))/(Isat**4.0))
print nSigma0Isat,Isat,cov
I0forFit=np.linspace(np.min(Raw2Avg),np.max(Raw2Avg),100)
data_fitted=isatOptThin(I0forFit,nSigma0Isat,Isat)

fig3=plt.figure()

pan3=fig3.add_subplot(1,1,1)
pan3.plot(Raw2Avg,Raw2Avg-Raw1Avg,'bo')
pan3.plot(I0forFit,data_fitted,'b-')
pan3.set_xlabel(r'$I_0$')
pan3.set_ylabel(r'$I_0-I_f$')
pan3.set_title(saveName +'\n'+r'$\sigma_0 n$=%.3f+/-%.3f, $I_{sat}$=%.0f +/- %.0f'%(nSigma0Isat/Isat,dnSigma0,Isat,dIsat))

fig3=plt.figure()
pan3=fig3.add_subplot(1,1,1)
pan3.plot(Raw2Avg,odAvg,'bo')
pan3.set_xlabel(r'$I_0$')
pan3.set_ylabel(r'od=-ln($I_f$/$I_0$)')