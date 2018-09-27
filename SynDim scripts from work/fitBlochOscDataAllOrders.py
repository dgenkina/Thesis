# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 13:53:37 2016

@author: dng5
"""

import Generic
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import numpy as np
import scipy.ndimage as snd
from scipy import optimize
import matplotlib.gridspec as gridspec
from sklearn.linear_model import Lasso, Ridge

def gaussLin3(x,amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C):
    gL = A+B*x+C*(x**2.0) + amp1*np.exp(-(x-x01)**2/sigma1**2) + amp2*np.exp(-(x-x02)**2/sigma2**2)+ amp3*np.exp(-(x-x03)**2/sigma3**2)
    return gL
    
def gaussLin3const(x,amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C):
    gL = A+B*x+C*(x**2.0) + amp1*np.exp(-(x-(x02-54.0))**2.0/sigma1**2.0) + amp2*np.exp(-(x-x02)**2.0/sigma2**2.0)+ amp3*np.exp(-(x-(x02+54.0))**2.0/sigma3**2.0)
    return gL
    
def plane((x,y), A,B,C):
    pl=A*x+B*y+C
    return pl.ravel()
    
def gauss2D((x,y),x0,y0,sigmaX,sigmaY,A,B):
    gauss=A*(np.exp(-(x-x0)**2.0/sigmaX**2.0-(y-y0)**2.0/sigmaY**2.0))+B
    return gauss.ravel()
    

def getRoi(array, image, xCent,yCent,r=16,eps=5,draw=True,color=(255,0,0)):
    ylen=array.shape[0]
    xlen=array.shape[1]
    bbox=(xCent-r,yCent-np.int(np.sqrt(r*r-eps*eps)),xCent+r,yCent+np.int(np.sqrt(r*r-eps*eps)))
    
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
    
    if draw: 
        draw = ImageDraw.Draw(image)
        draw.ellipse(bbox,outline=color)
    return (counts, countsL, countsR, countsT, countsB)
    


roi=np.array([360, 670, 410,630])
roiFlea=([280,330,200,250])#([205,260,365,400])
filestart=4
filestop=164
fileroot = 'Y:/Data/2016/August/24/PIXIS_24Aug2016'  
filerootFlea='X:/2016/August/24/Flea3_24Aug2016'  
filelist=np.arange(124,184)#np.append(np.arange(154,182),np.arange(186,204))#
key='pulseDelay'
date=fileroot.split("/")[-1].split("_")[1]
saveName=date+'_files_'+np.str(filelist[0])+'-'+str(filelist[-1])+'_allOrders'
print 'one'

(fractionP,fraction0,fractionM,waveDict,imbalArray,signalGood, qlist, xCenters, fracAll)=blochOscFractionsV3(fileroot,filelist,roi,key,plot=True,xlabel='Oscillation time [s]',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea,weight=False)        
np.savez(saveName,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,imbalArray=imbalArray,qlist=qlist,signalGood=signalGood,tlist=waveDict[key],xCenters=xCenters,fracAll=fracAll)
#

def fringeremoval(filelist,filenum, roi, fileroot=fileroot):
    angle=-41
    filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=angle)   
    od=snd.interpolation.rotate(dict1['OptDepth'],angle)[roi[0]:roi[1],roi[2]:roi[3]]   
    Raw1=snd.interpolation.rotate(dict1['Raw1'],angle)[roi[0]:roi[1],roi[2]:roi[3]]
    A=np.array(Raw1.flatten())
    R=[]
    for num in filelist:
        filename=fileroot+"_"+ str(num).zfill(4) + ".ibw"
        dict1 =processIBW(filename, angle=-41)  
        R.append(snd.interpolation.rotate(dict1['Raw2'],angle)[roi[0]:roi[1],roi[2]:roi[3]].flatten())
    R=np.array(R).transpose()
    print R.shape
    print A.shape   
#        
#    lasso = Lasso(alpha=0.01,max_iter=100000)
#    lasso.fit(R, A)
#    c = lasso.coef_

    ridge = Ridge(alpha=0.01)
    ridge.fit(R, A)
    c = ridge.coef_
    
    print c
    Raw2fitted=np.dot(R,c).reshape(Raw1.shape)
    odClean = -np.log(Raw1.astype(float)/Raw2fitted.astype(float))
    odClean = readIgor.zeroNansInfsVector(odClean)
#    
#    figure=plt.figure()
#    pan=figure.add_subplot(1,2,1)
#    pan.imshow(od, vmin=-0.15, vmax=0.3)
#    pan2=figure.add_subplot(1,2,2)
#    pan2.imshow(odClean, vmin=-0.15, vmax=0.3)
    return odClean

#cutoff=0.05
#fieldGoodArray=((imbalArray<cutoff) & signalGo)
#fractionP=fractionP[fieldGoodArray]
#fraction0=fraction0[fieldGoodArray]
#fractionM=fractionM[fieldGoodArray]
#time=waveDict[key][fieldGoodArray]
#figure=plt.figure()
#panel=figure.add_subplot(1,1,1)
#panel.plot(time,fractionP,'bo', label=r'$m_F$=+1')
#panel.plot(time,fraction0,'go', label=r'$m_F$=0')
#panel.plot(time,fractionM,'ro', label=r'$m_F$=-1')
#panel.set_xlabel(key)
#legend()

def checkFleaRoi(filerootFlea,filenum,roiFlea):
    filenameFlea=filerootFlea+"_"+ str(filenum).zfill(4) + ".ibw"
    dictFlea =processIBW(filenameFlea, angle=-38)
    od1=dictFlea['od1'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]
    od2=dictFlea['od2'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]   
    fig=plt.figure()
    pan1=fig.add_subplot(2,2,1)
    pan1.imshow(od1)
    pan2=fig.add_subplot(2,2,2)
    pan2.imshow(od2)
    pan3=fig.add_subplot(2,2,3)
    pan3.imshow(dictFlea['od1'])
    pan4=fig.add_subplot(2,2,4)
    pan4.imshow(dictFlea['od2'])
    return
    
def blochOscOneFileV3(fileroot,filenum, roi, draw=True, checkField=False,filerootFlea=fileroot,roiFlea=roi,bgndRidge=False,weight=True):
    filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=-44)
    print filename
    
    if checkField:
        roiB=np.array([roiFlea[0],roiFlea[1],2*roiFlea[2]-roiFlea[3],roiFlea[2]])
        filenameFlea=filerootFlea+"_"+ str(filenum).zfill(4) + ".ibw"
        dictFlea =processIBW(filenameFlea, angle=-41)
        od1=dictFlea['od1'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]
        od1B=dictFlea['od1'][roiB[0]:roiB[1],roiB[2]:roiB[3]]
        od2=dictFlea['od2'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]   
        od2B=dictFlea['od2'][roiB[0]:roiB[1],roiB[2]:roiB[3]]
        num1=np.sum(od1)-np.sum(od1B)
        num2=np.sum(od2)-np.sum(od2B)
        imbal=(num1-num2)/(num1+num2)
        signalGood=((num1>0) & (num2>0))
        print num1, num2, imbal, signalGood
    else:
        imbal = nan
        signalGood=False
    
    
    if bgndRidge:
        odRoi=fringeremoval(filelist,filenum, roi, fileroot=fileroot)
        
    else:
        odRoi=dict1['rotODcorr'][roi[0]:roi[1],roi[2]:roi[3]]   
    
    delta=30
    xM=np.array([])
    yM=np.array([])
    fM=np.array([])
    xvals = np.arange(roi[2]-delta,roi[3]+delta)
    yvals = np.arange(roi[0]-delta,roi[1]+delta)
    for x in xvals:
        for y in yvals:
            if (x>roi[2] and x<roi[3] and y>roi[0] and y<roi[1])==False:
                xM=np.append(xM,x)
                yM=np.append(yM,y)
                fM=np.append(fM,dict1['rotODcorr'][y,x])
                
    (A,B,C), cov = optimize.curve_fit(plane,(xM,yM),fM,p0=(0.03,(roi[3]-roi[2])/2.0,(roi[1]-roi[0]/2.0)))
#    print A,B,C,D,F,G
    (x,y)=np.meshgrid(np.arange(roi[2],roi[3]),np.arange(roi[0],roi[1]))
    fitted_plane=plane((x,y),A,B,C).reshape(np.arange(roi[0],roi[1]).size,np.arange(roi[2],roi[3]).size)
    odRoi1=odRoi-fitted_plane
    
    odFiltered=snd.filters.gaussian_filter(odRoi1,1.0)
    
    if draw:      
        fig1=plt.figure()
        pan1=fig1.add_subplot(1,1,1)
        pan1.imshow(odFiltered,vmin=-0.15,vmax=0.3)
        
    
    odFiltLin = np.sum(odFiltered,axis=0)
        
    peak2=np.max(odFiltLin[100:150])

    xGuess2=np.float( np.where(odFiltLin==peak2)[0])
#    xGuess2=175
#    xGuess1=xGuess2-55.0
#    xGuess3=xGuess2+52.0
#    
    if draw:
        figure2=plt.figure()
        panel2=figure2.add_subplot(1,1,1)
        panel2.plot(odFiltLin,'bo')
    
#    p=(1.5,xGuess1,10.0,peak2,xGuess2,10.0,0.5,xGuess3,10.0,0,0,0)
#    (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3, np.arange(odFiltLin.size), odFiltLin, p0 =p )
#    data_fitted = gaussLin3(np.arange(odFiltLin.size), *(amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C))
    
#    if draw:
#        panel2.plot(data_fitted,'b-')
  
    xCent=xGuess2

    
    roi2=np.array([620, 880, 630, 870])
    

    
    rotOD2=snd.interpolation.rotate(dict1['rotOD'],50)
    odRoiRot2 = rotOD2[roi2[0]:roi2[1],roi2[2]:roi2[3]]

    
    delta=20
    xM=np.array([])
    yM=np.array([])
    fM=np.array([])
    xvals = np.arange(roi2[2]-delta,roi2[3]+delta)
    yvals = np.arange(roi2[0]-delta,roi2[1]+delta)
    for x in xvals:
        for y in yvals:
            if (x>roi2[2] and x<roi2[3] and y>roi2[0] and y<roi2[1])==False:
                xM=np.append(xM,x)
                yM=np.append(yM,y)
                fM=np.append(fM,rotOD2[y,x])
                
    (A,B,C), cov = optimize.curve_fit(plane,(xM,yM),fM)
    
    (x,y)=np.meshgrid(np.arange(roi2[2],roi2[3]),np.arange(roi2[0],roi2[1]))
    fitted_plane2=plane((x,y),A,B,C).reshape(np.arange(roi2[0],roi2[1]).size,np.arange(roi2[2],roi2[3]).size)
    odRoi2=odRoiRot2-fitted_plane2
    
    odFiltered2=snd.filters.gaussian_filter(odRoi2,1.0)
    
    if draw:
        fig2=plt.figure()
        pan2=fig2.add_subplot(1,1,1)
        pan2.imshow(odFiltered2 , vmin=-0.15, vmax=0.3)
        
        
    odFiltLin2 = np.sum(odFiltered2,axis=0)   

    peak2=np.max(odFiltLin2[100:150])
    xGuess2=np.float(np.where(odFiltLin2==peak2)[0])


#    xGuess1=xGuess2-55.0
#    xGuess3=xGuess2+54.0
    
    if draw:     
        figure3=plt.figure()
        panel3=figure3.add_subplot(1,1,1)
        panel3.plot(odFiltLin2,'bo')
        
#    p=(1.0,10.0,1.0,xGuess2,10.0,1.0,10.0,4.0,0,0)
#    (amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3const, np.arange(odFiltLin2.size), odFiltLin2,p0=p )
#    data_fitted = gaussLin3const(np.arange(odFiltLin2.size), *(amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C))

    
   # yCent=0.865*xGuess2+54 #odtkick negative
  #  yCent=1.32*xGuess2-3.3# odtkick positive
  #  yCent=-0.055*xGuess2**2.0+13.78*xGuess2-690.372

#    print xCent
#    print xGuess2
    
    

    yCent=-0.92*xCent+1.4*xGuess2+120.0
    
        
#    print yCent
#    print xCent
#    xCent=165
    
    print xCent,xGuess2, yCent
#    print yCent
    #yCent=118
    #xCent=100
    
    norm=mpl.colors.Normalize(vmin=-0.15,vmax=0.3)
    im = Image.fromarray(np.uint8(plt.cm.jet(norm(odFiltered))*255))
    y1=71
    x1=59
    y2=23
    x2=2
    yrf=2
    offsets = np.array([[0,0],[-x2,y1],[x2,-y1],[-x2*2,y1*2],[x2*2,-y1*2],[-x1,-y2],[-x1-x2,-y2+y1],[-x1+x2,-y2-y1],[-x1-x2*2,-y2+y1*2],[-x1+x2*2,-y2-y1*2],[x1,y2],[x1-x2,y2+y1],[x1+x2,y2-y1],[x1-x2*2,y2+y1*2],[x1+x2*2,y2-y1*2]])#np.array([[0,0],[0,69],[0,-69],[-58,49],[-58,119],[-58,-21],[57,-47],[57,21],[57,-115]])
    offsetsShort=   np.array([[0,0],[0,y1],[0,-y1],[-x1,-y2],[-x1,-y2+y1],[-x1,-y2-y1],[x1,y2],[x1,y2+y1],[x1,y2-y1]]) 
 
  #  offsets=offsetsShort
    counts=np.zeros(offsets.shape[0])    
    for i in np.arange(offsets.shape[0]):
        offs=offsets[i]
        counts[i]= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=17,draw=False)[0]
#    print counts
    maxInd=np.where(counts==np.max(counts))[0][0]
    
    allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],r=17)
#    print allcounts
    i=0
    while ((np.abs(allcounts[4]-allcounts[3])>2.0) & (i<20)):
        if (allcounts[4]-allcounts[3])>0:
            yCent=yCent+1
       #     print "new yCent = " +str(yCent)
            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
        else: 
            yCent=yCent-1
      #      print "new yCent = " +str(yCent)
            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
        i=i+1
#    print i
    i=0
    while ((np.abs(allcounts[2]-allcounts[1])>2.0) & (i<20)):
        if (allcounts[2]-allcounts[1])>0:
            xCent=xCent+1
     #       print "new xCent = " +str(xCent)blochOscFractionsV2(fileroot,filelist,roi,key,plot=True,xlabel='',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea
            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
        else: 
            xCent=xCent-1
     #       print "new xCent = " +str(xCent)a[6]['Note']
            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
        i=i+1
#    print i
 #   print allcounts
    
    for i in np.arange(offsets.shape[0]):
        offs=offsets[i]
        r=9.0    
        if int(i/(offsets.size/3))==0:
            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
#            if count>0:
#                count -= 6.0
          #  bgnd1=getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1]+2.0*r,r=r, draw=draw,color=(0,0,0))[0]
            bgnd=getRoi(odFiltered, im, xCent+offs[0]-2.0*r,yCent+offs[1],r=r, draw=draw,color=(0,0,0))[0]
            counts[i]=count-bgnd
            if weight:
                counts[i]=counts[i]/((1.6e-5)*offs[1]**2.0-0.00583*offs[1]+1.296)
        if int(i/(offsets.size/3))==1:
            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
#            if count>0:
#                count -= 6.0
            bgnd=getRoi(odFiltered, im, xCent+offs[0]-2.0*r,yCent+offs[1],r=r, draw=draw,color=(0,0,0))[0]
            counts[i]=count-bgnd
            if weight:
                counts[i]=counts[i]/((1.5e-5)*offs[1]**2.0-0.00538*offs[1]+1.329)
        else:
            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
#            if count>0:
#                count -= 6.0
            bgnd=getRoi(odFiltered, im, xCent+offs[0]-2.0*r,yCent+offs[1],r=r, draw=draw,color=(0,0,0))[0]
            counts[i]=count-bgnd
            if weight:
                counts[i]=counts[i]/((1.4e-5)*offs[1]**2.0-0.00473*offs[1]+1.266)
    if draw:
        fig4=plt.figure()
        pan4=fig4.add_subplot(1,1,1)
        pan4.imshow(im)
        pan4.set_title(filename)
        
    return counts, odFiltered, xCent,yCent, imbal,signalGood,  dict1
    
def blochOscFractionsV3(fileroot,filelist,roi,key,plot=True,xlabel='',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea,weight=True):
   # filerange = np.arange(filestart,filestop+1)    

    fractionP=np.zeros(filelist.size)
    fraction0=np.zeros(filelist.size)
    fractionM=np.zeros(filelist.size)
    imbalArray=np.ones(filelist.size)
    qlist=np.zeros(filelist.size)
    signalGood=np.ones(filelist.size,dtype=bool)
    
#    unindexedWave=np.zeros(filelist.size)
   # filename = fileroot + "_"+ str(filelist[0]).zfill(4) + ".ibw"

    (counts, odFiltered, xCent,yCent,imbal,signalG, dict1)=blochOscOneFileV2(fileroot,filelist[0], roi, draw=False)
    fracAll=np.zeros((filelist.size,counts.size))
    infoString=dict1["Note"]
    waveDict=getIndexedWaves(infoString)
    waveDict=waveDict.fromkeys(waveDict.keys(),[])
    odFiltAll=np.zeros((filelist.size,odFiltered.shape[0],odFiltered.shape[1]))
    xCenters=np.zeros(filelist.size)
    for ind, filenum in enumerate(filelist):
       # filename = fileroot + "_"+ str(filenum).zfill(4) + ".ibw"
        #print filename
        (counts, odFiltered, xCent, yCent,imbal, signalG,dict1)=blochOscOneFileV2(fileroot,filenum, roi, draw=False,checkField=checkField,filerootFlea=filerootFlea,roiFlea=roiFlea,weight=weight)
        print filenum
        xCenters[ind]=xCent
        qlist[ind]=2.0*(103.5-yCent)/71.0
        if checkField:
            imbalArray[ind]=imbal
            signalGood[ind]=signalG
        infoString=dict1["Note"]
        waveDictLocal=getIndexedWaves(infoString)
        for k in waveDict.iterkeys():
            
            waveDict[k]=np.append(waveDict[k],waveDictLocal[k])

 #       unindexedWave[ind]=Generic.StringByKey('nan', key, infoString, '=', '\n')
        odFiltAll[ind]=odFiltered
        
        fracAll[ind]=counts/np.sum(counts)
        orders=counts.size/3
        fractionP[ind]=np.sum(fracAll[ind][orders:2*orders])
        fraction0[ind]=np.sum(fracAll[ind][:orders])
        fractionM[ind]=np.sum(fracAll[ind][2*orders:])
    if plot:
        figure3=plt.figure()
        panel3=figure3.add_subplot(1,1,1)
        if key == 'ind':
            xvals=range(filelist.size)
        else:
            xvals =waveDict[key]
            
        if xlabel=='':
            if key == 'ind':
                xlabel = 'shot number'
            else:
                xlabel = key
        panel3.plot(xvals,fractionP,'bo', label='mF=+1')
        panel3.plot(xvals,fraction0,'go', label='mF=0')
        panel3.plot(xvals,fractionM,'ro', label='mF=-1')
        panel3.set_xlabel(xlabel)
        panel3.set_ylabel('Fractional population')
        legend()
        
        figure4=plt.figure()
        pan4=figure4.add_subplot(1,1,1)
        pan4.plot(xvals, qlist,'bo')
        pan4.set_ylabel(r'Quasimomentum [$k_L$]')
        pan4.set_xlabel('Oscillation time [s]')
        
        figure = plt.figure()
        gs = gridspec.GridSpec(3,filelist.size)
        gs.update(left=0.01, right=0.99, wspace=0.005, hspace=0)
        if key == 'ind':
            sort=range(filelist.size)
            tsorted=range(filelist.size)
        else:
            sort=np.argsort(waveDict[key])
            tsorted = waveDict[key][sort]

        w=30
        for i in range(filelist.size):
            j=sort[i]
            panel1=figure.add_subplot(gs[i])
            panel1.imshow(odFiltAll[j][:,int(xCenters[j]-49-w):int(xCenters[j]-49+w)],vmin=-0.05, vmax=0.3)
            panel1.yaxis.set_ticks([])  
            panel1.xaxis.set_ticks([])
    
            panel2=figure.add_subplot(gs[i+filelist.size])
            panel2.imshow(odFiltAll[j][:,int(xCenters[j]-w):int(xCenters[j]+w)],vmin=-0.05, vmax=0.3)
            panel2.yaxis.set_ticks([])  
            panel2.xaxis.set_ticks([])
            
            panel3=figure.add_subplot(gs[i+2*filelist.size])
            panel3.imshow(odFiltAll[j][:,int(xCenters[j]+49-w):int(xCenters[j]+49+w)],vmin=-0.05, vmax=0.3)
            panel3.yaxis.set_ticks([])  
            panel3.xaxis.set_ticks([])      
        
    return fractionP,fraction0,fractionM,waveDict,imbalArray,signalGood, qlist, xCenters, fracAll