# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 17:27:52 2016

@author: dng5
"""

#import readIgor, Generic
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import numpy as np
import scipy.ndimage as snd
from scipy import optimize
import matplotlib.gridspec as gridspec
from sklearn.linear_model import Lasso, Ridge
import time


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
    


def getRoi(array, image, xCent,yCent,weightArray=1.0,updateWeights=False,r=16,eps=5,draw=True,color=(255,0,0),horizontalStretch=True):
    ylen=array.shape[0]
    xlen=array.shape[1]
    if horizontalStretch:
        bbox=(xCent-r,yCent-np.int(np.sqrt(r*r-eps*eps)),xCent+r,yCent+np.int(np.sqrt(r*r-eps*eps)))
    else:
        bbox=(xCent-np.int(np.sqrt(r*r-eps*eps)),yCent-r,xCent+np.int(np.sqrt(r*r-eps*eps)),yCent+r)
    
    y,x = np.ogrid[-yCent:ylen-yCent, -xCent:xlen-xCent]
    if horizontalStretch:
        xloc1=x-eps
        yloc1=y
        xloc2=x+eps
        yloc2=y
    else:
        xloc1=x
        yloc1=y-eps
        xloc2=x
        yloc2=y+eps
    
    mask = np.sqrt((xloc1)**2.0 + yloc1**2.0) +np.sqrt((xloc2)**2.0 + yloc2**2.0) <= 2.0*r
    if updateWeights:
        weightArray[mask]=0
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
    return (counts, countsL, countsR, countsT, countsB, weightArray)
    

roi=np.array([380, 690, 350,700])
roiFlea=([280,330,200,250])#([205,260,365,400])
filestart=4
filestop=164
fileroot = 'Y:/Data/2017/March/09/PIXIS_09Mar2017'  
filerootFlea= 'Y:/Data/2017/March/09/Flea3_09Mar2017' 
filelist=np.arange(754,784)#np.append(np.arange(634,637),np.arange(639,664))#
key='pulseDelay'
date=fileroot.split("/")[-1].split("_")[1]
saveName=date+'_files_'+np.str(filelist[0])+'-'+str(filelist[-1])#+'_weighted'
checkCoherence=True
w0=0.8
wP=1.12
wM=1.0
filelistProbes=np.arange(694, 794)
fitProbes=True

def getProbeReconstruction(fileroot, filelistProbes, roi):
    filename=fileroot+"_"+ str(filelistProbes[0]).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=-43)
    rotProbe=snd.interpolation.rotate(dict1['Raw2'],-43)[roi[0]:roi[1],roi[2]:roi[3]]
    N=rotProbe.flatten().size
    print 'number of pixels in roi = ' + str(N)
    t1=time.clock()
    probeMatrix=np.ones((N,filelistProbes.size+1))
    for ind, filenum in enumerate(filelistProbes):
        filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
        dict1 =processIBW(filename, angle=-43)
        rotProbe=dict1['rotRaw2'][roi[0]:roi[1],roi[2]:roi[3]]
        probeMatrix[:,ind+1]=rotProbe.flatten()
    probeMatrixT=probeMatrix.transpose()
    t2=time.clock()
    print 'constructed probe matrix in'+str(t2-t1)
    return probeMatrix, probeMatrixT
#    
#if fitProbes:
#    probeMatrix,probeMatrixT=getProbeReconstruction(fileroot, filelistProbes, roi)
#else:
#    probeMatrix=1.0
#    probeMatrixT=1.0

(numM, num0, numP, fractionM, fraction0, fractionP, waveDict,num1Array,num2Array,imbalArray,signalGood,qlist,xCenters,coherentList)=blochOscFractionsV2Rf(fileroot,filelist,roi,key,plot=True,xlabel='Oscillation time [s]',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea,weights=False)        
np.savez(saveName,num0=num0,numM=numM,numP=numP,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,num1Array=num1Array,num2Array=num2Array,imbalArray=imbalArray,xCenters=xCenters,signalGood=signalGood,tlist=waveDict[key],qlist=qlist,coherentList=coherentList)
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
    odClean = zeroNansInfsVector(odClean)
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
IsatCounts = 477.6
imagingTime=40.0

def zeroNansInfs(x):
    if x!=x:
        return 0
    elif np.isinf(x):
        return 0
    else:
        return x
zeroNansInfsVector = np.vectorize(zeroNansInfs, otypes=[np.float])

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
def blochOscOneFileV2Rf(fileroot,filenum, roi, draw=True, checkField=False,filerootFlea=fileroot,roiFlea=roi,bgndRidge=False):
    filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=-43)
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
        num1=nan
        num2=nan
    
    
    if bgndRidge:
        odRoi=fringeremoval(filelist,filenum, roi, fileroot=fileroot)
        
    else:
        odRoi=dict1['rotODcorr'][roi[0]:roi[1],roi[2]:roi[3]]   
        atomsRoi=dict1['rotRaw1'][roi[0]:roi[1],roi[2]:roi[3]]  
    
#    delta=30
#    xM=np.array([])
#    yM=np.array([])
#    fM=np.array([])
#    xvals = np.arange(roi[2]-delta,roi[3]+delta)
#    yvals = np.arange(roi[0]-delta,roi[1]+delta)
#    for x in xvals:
#        for y in yvals:
#            if (x>roi[2] and x<roi[3] and y>roi[0] and y<roi[1])==False:
#                xM=np.append(xM,x)
#                yM=np.append(yM,y)
#                fM=np.append(fM,dict1['rotODcorr'][y,x])
#                
    (x,y)=np.meshgrid(np.arange(roi[2],roi[3]),np.arange(roi[0],roi[1]))
    
    (A,B,C), cov = optimize.curve_fit(plane,(x,y),odRoi.flatten(),p0=(0.0,0.0,0.0))
    print A,B,C    
    fitted_plane=plane((x,y),A,B,C).reshape(np.arange(roi[0],roi[1]).size,np.arange(roi[2],roi[3]).size)
    odRoi1=odRoi-fitted_plane
    
    odFiltered=snd.filters.gaussian_filter(odRoi1,1.0)#+buildOD
    
    if draw:      
        fig1=plt.figure()
        pan1=fig1.add_subplot(1,1,1)
        pan1.imshow(odFiltered,vmin=-0.15,vmax=0.3)
        
    
    odFiltLin = np.sum(odFiltered,axis=0)
        
#    peakL=np.max(odFiltLin[:120])
#    xGuessL=np.float( np.where(odFiltLin==peakL)[0])
#    print xGuessL
#    peakR=np.max(odFiltLin[120:])
#
#    xGuessR=np.float( np.where(odFiltLin==peakR)[0])
#    print xGuessR
#    peakT=np.max(odFiltLin)blochOscOneFileV2Rf(fileroot,filenum, roi
#    xGuessT=np.float( np.where(odFiltLin==peakT)[0])
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
  
#    xCent=(xGuessR+xGuessL)/2.0

    peak=np.max(odFiltLin[140:180])
    xGuess=np.float(np.where(odFiltLin==peak)[0])
    xCent=xGuess    
    


#    xGuess1=xGuess2-55.0
#    xGuess3=xGuess2+54.0

        
#    p=(1.0,10.0,1.0,xGuess2,10.0,1.0,10.0,4.0,0,0)
#    (amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3const, np.arange(odFiltLin2.size), odFiltLin2,p0=p )
#    data_fitted = gaussLin3const(np.arange(odFiltLin2.size), *(amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C))

    
   # yCent=0.865*xGuess2+54 #odtkick negative
  #  yCent=1.32*xGuess2-3.3# odtkick positive
  #  yCent=-0.055*xGuess2**2.0+13.78*xGuess2-690.372


    
        
    #for rf only:
    odFiltVert = np.sum(odFiltered,axis=1)  
    plot(odFiltVert)
    peakRf=np.max(odFiltVert)
    yCent = np.float( np.where(odFiltVert==peakRf)[0])
    print yCent
    print xCent
#    xCent=165
#    yCent=183
  #  print xCent,xGuess2, yCent
#    print yCent
#    yCent=135
    
    norm=mpl.colors.Normalize(vmin=-0.15,vmax=0.3)
    im = Image.fromarray(np.uint8(plt.cm.jet(norm(odFiltered))*255))
    y1=61
    x1=50
    yrf=0
 #   offsets = np.array([[0,0],[0,y1],[0,-y1],[0,y1*2],[0,-y1*2],[-x1,-y2],[-x1,-y2+y1],[-x1,-y2-y1],[-x1,-y2+y1*2],[-x1,-y2-y1*2],[x1,y2],[x1,y2+y1],[x1,y2-y1],[x1,y2+y1*2],[x1,y2-y1*2]])#np.array([[0,0],[0,69],[0,-69],[-58,49],[-58,119],[-58,-21],[57,-47],[57,21],[57,-115]])
 #   offsetsShort=   np.array([[0,0],[0,y1],[0,-y1],[-x1,-y2],[-x1,-y2+y1],[-x1,-y2-y1],[x1,y2],[x1,y2+y1],[x1,y2-y1]]) 
    rfOffsets = np.array([[0,0],[0,y1],[0,-y1],[-x1,yrf],[-x1,y1+yrf],[-x1,-y1+yrf],[x1,-yrf],[x1,y1-yrf],[x1,-y1-yrf]])    
    rfOffsetsShortest=np.array([[0,0],[-x1,yrf],[x1,-yrf]])
    offsets=rfOffsets
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
    r=16.0 
    if fitProbes:
        weightArray=np.ones(odFiltered.shape)
    for i in np.arange(offsets.shape[0]):
        offs=offsets[i]

        count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
        if fitProbes:
            (count,cL,cR,cT,cB,weightArray)= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r,weightArray=weightArray,updateWeights=True, draw=draw)
            

        bgnd1=getRoi(odFiltered, im, xCent+offs[0]+2.0*r,yCent+offs[1],r=r, draw=draw,color=(0,0,0))[0]
        bgnd2=getRoi(odFiltered, im, xCent+offs[0]-2.0*r,yCent+offs[1],r=r, draw=draw,color=(0,0,0))[0]
        counts[i]=count-(bgnd1+bgnd2)/2.0
        
    if fitProbes:
        xTw=probeMatrixT*weightArray.flatten()
        rhs=np.dot(xTw,atomsRoi.flatten())

        lhs=np.dot(xTw,probeMatrix)

        beta=np.linalg.solve(lhs,rhs)
        newProbe=np.dot(probeMatrix,beta).reshape(atomsRoi.shape)
        newOd=-np.log(atomsRoi/newProbe)
        newOd = zeroNansInfsVector(newOd)
        odFiltered = newOd+(newProbe-atomsRoi)/(IsatCounts*imagingTime)
        odFiltered=snd.filters.gaussian_filter(odRoi1,1.0)
        
        im2 = Image.fromarray(np.uint8(plt.cm.jet(norm(odFiltered))*255))
        for i in np.arange(offsets.shape[0]):
            offs=offsets[i]
            count= getRoi(odFiltered, im2, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
            counts[i]=count
        


    if checkCoherence:
        maxInd=np.int(maxInd/3)*3
        
        extraAtoms1=getRoi(odFiltered,im,xCent+(offsets[maxInd,0]+offsets[maxInd+1,0])/2.0,yCent+(offsets[maxInd,1]+offsets[maxInd+1,1])/2.0,r=r)[0]#r=y1/2.0-r,eps=np.sqrt((y1/2.0-r)**2.0-r**2.0),horizontalStretch=False)[0]
        extraAtoms1bgnd=getRoi(odFiltered,im,xCent+(offsets[maxInd,0]+offsets[maxInd+1,0])/2.0+2.0*r,yCent+(offsets[maxInd,1]+offsets[maxInd+1,1])/2.0,r=r,color=(0,0,0))[0]#r=y1/2.0-r,eps=np.sqrt((y1/2.0-r)**2.0-r**2.0),color=(0,0,0),horizontalStretch=False)[0]
        extraAtoms2=getRoi(odFiltered,im,xCent+(offsets[maxInd,0]+offsets[maxInd+2,0])/2.0,yCent+(offsets[maxInd,1]+offsets[maxInd+2,1])/2.0,r=r)[0]#r=y1/2.0-r,eps=np.sqrt((y1/2.0-r)**2.0-r**2.0),horizontalStretch=False)[0]
        extraAtoms2bgnd=getRoi(odFiltered,im,xCent+(offsets[maxInd,0]+offsets[maxInd+2,0])/2.0+2.0*r,yCent+(offsets[maxInd,1]+offsets[maxInd+2,1])/2.0,r=r,color=(0,0,0))[0]#r=y1/2.0-r,eps=np.sqrt((y1/2.0-r)**2.0-r**2.0),color=(0,0,0),horizontalStretch=False)[0]
        extraAtoms1=extraAtoms1-extraAtoms1bgnd
        extraAtoms2=extraAtoms2-extraAtoms2bgnd
        print 'Extra atoms:'
        print extraAtoms1,extraAtoms2
        print counts[maxInd], counts[maxInd+1], counts [maxInd+2]
        goodAtoms=counts[maxInd]+ counts[maxInd+1]+counts [maxInd+2]
        extraAtoms=np.max([extraAtoms1,extraAtoms2])
        if extraAtoms>goodAtoms/3.0:
            coherent=False
        elif goodAtoms<0.0:
            coherent=False
        else:
            coherent=True
    else:
        coherent=True

    if draw:
        fig4=plt.figure()
        pan4=fig4.add_subplot(1,1,1)
        pan4.imshow(im)
        pan4.set_title(filename)
        if fitProbes:
            fig4=plt.figure()
            pan4=fig4.add_subplot(1,1,1)
            pan4.imshow(im2)
            pan4.set_title(filename+'_reconProbe')
    
    print 'coherent = '
    print coherent
    
        
    return counts, odFiltered, xCent,yCent,num1,num2, imbal,signalGood,coherent,  dict1
    
def blochOscFractionsV2Rf(fileroot,filelist,roi,key,plot=True,xlabel='',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea,weights=False):
   # filerange = np.arange(filestart,filestop+1)    

    num0=np.zeros(filelist.size)
    numM=np.zeros(filelist.size)
    numP=np.zeros(filelist.size)

    fractionP=np.zeros(filelist.size)
    fraction0=np.zeros(filelist.size)
    fractionM=np.zeros(filelist.size)
    
    num1Array=np.ones(filelist.size)
    num2Array=np.ones(filelist.size)
    imbalArray=np.ones(filelist.size)
    qlist=np.zeros(filelist.size)
    signalGood=np.ones(filelist.size,dtype=bool)
    coherentList=np.ones(filelist.size,dtype=bool)
   # unindexedWave=np.zeros(filelist.size)
   # filename = fileroot + "_"+ str(filelist[0]).zfill(4) + ".ibw"

    (counts, odFiltered, xCent,yCent,num1,num2,imbal,signalGo, coherent,dict1)=blochOscOneFileV2Rf(fileroot,filelist[0], roi,draw=False)
    infoString=dict1["Note"]
    waveDict=getIndexedWaves(infoString)
    waveDict=waveDict.fromkeys(waveDict.keys(),[])
    odFiltAll=np.zeros((filelist.size,odFiltered.shape[0],odFiltered.shape[1]))
    xCenters=np.zeros(filelist.size)
    yCenters=np.zeros(filelist.size)
    for ind, filenum in enumerate(filelist):
       # filename = fileroot + "_"+ str(filenum).zfill(4) + ".ibw"
        #print filename
        (counts, odFiltered, xCent, yCent,num1,num2,imbal, signalG,coherent, dict1)=blochOscOneFileV2Rf(fileroot,filenum, roi,draw=False,checkField=checkField,filerootFlea=filerootFlea,roiFlea=roiFlea)
        print filenum
        xCenters[ind]=xCent
        yCenters[ind]=yCent
        qlist[ind]=2.0*(103.5-yCent)/61.0
        if checkField:
            imbalArray[ind]=imbal
            signalGood[ind]=signalG
            coherentList[ind]=coherent
            num1Array[ind]=num1
            num2Array[ind]=num2
        infoString=dict1["Note"]
        waveDictLocal=getIndexedWaves(infoString)
        for k in waveDict.iterkeys():
            
            waveDict[k]=np.append(waveDict[k],waveDictLocal[k])

  #      unindexedWave[ind]=Generic.StringByKey('nan', key, infoString, '=', '\n')
        odFiltAll[ind]=odFiltered
        
#
#        roiNum=counts.size/3
#        num1[ind]=np.sum(counts[roiNum:roiNum*2])
#        if num1[ind]<0.0:
#            print 'num1 < 0 in file ' + str(filenum)
#        num2[ind]=np.sum(counts[0:roiNum])
#        if num2[ind]<0.0:
#            print 'num2 < 0 in file ' + str(filenum)
#        num3[ind]=np.sum(counts[roiNum*2:roiNum*3])
#        if num3[ind]<0.0:
#            print 'num3 < 0 in file ' + str(filenum)
        roiNum=counts.size/3
        num0[ind]=np.sum(counts[0:roiNum])
        numP[ind]=np.sum(counts[roiNum:2*roiNum])
        numM[ind]=np.sum(counts[roiNum*2:roiNum*3])

        total=num0[ind]+numM[ind]+numP[ind]
        
        if weights:
            num0[ind]=w0*num0[ind]
            numP[ind]=wP*numP[ind]
            numM[ind]=wM*numM[ind]   
            total=num0[ind]+numM[ind]+numP[ind]
        

        fractionP[ind]=numP[ind]/total
        fraction0[ind]=num0[ind]/total
        fractionM[ind]=numM[ind]/total

        print numM[ind], num0[ind], numP[ind]
        print fractionM[ind], fraction0[ind], fractionP[ind]
        
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

        panel3.set_xlabel(key)
        panel3.set_ylabel('Fractional population')
        legend()
        
        figure4=plt.figure()
        pan4=figure4.add_subplot(1,1,1)
        pan4.plot(xvals, qlist,'bo')
        pan4.set_ylabel(r'Quasimomentum [$k_L$]')
        pan4.set_xlabel(key)
        
        figure = plt.figure()
        gs = gridspec.GridSpec(3,filelist.size)
        gs.update(left=0.01, right=0.99, wspace=0.005, hspace=0)
        if key == 'ind':
            sort=range(filelist.size)
            tsorted=range(filelist.size)
            waveDict['ind']=tsorted
        else:
            sort=np.argsort(waveDict[key])
            tsorted = waveDict[key][sort]
         
        w=30
        for i in range(filelist.size):
            j=sort[i]
            panel2=figure.add_subplot(gs[i])
            panel2.imshow(odFiltAll[j][:,int(xCenters[j]-48-w):int(xCenters[j]-48+w)],vmin=-0.05, vmax=0.3)
            panel2.yaxis.set_ticks([])  
            panel2.xaxis.set_ticks([])
    
            panel3=figure.add_subplot(gs[i+filelist.size])
            panel3.imshow(odFiltAll[j][:,int(xCenters[j]-w):int(xCenters[j]+w)],vmin=-0.05, vmax=0.3)
            panel3.yaxis.set_ticks([])  
            panel3.xaxis.set_ticks([])
            
            panel4=figure.add_subplot(gs[i+2*filelist.size])
            panel4.imshow(odFiltAll[j][:,int(xCenters[j]+48-w):int(xCenters[j]+48+w)],vmin=-0.05, vmax=0.3)
            panel4.yaxis.set_ticks([])  
            panel4.xaxis.set_ticks([])  
 
        
    return numM, num0, numP, fractionM, fraction0, fractionP, waveDict,num1Array,num2Array,imbalArray,signalGood,qlist,xCenters,coherentList    
#def blochOscOneFileV2(fileroot,filenum, roi, draw=True, checkField=False,filerootFlea=fileroot,roiFlea=roi,bgndRidge=False):
#    filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
#    dict1 =processIBW(filename, angle=-44)
#    print filename
#    
#    if checkField:
#        roiB=np.array([roiFlea[0],roiFlea[1],2*roiFlea[2]-roiFlea[3],roiFlea[2]])
#        filenameFlea=filerootFlea+"_"+ str(filenum).zfill(4) + ".ibw"
#        dictFlea =processIBW(filenameFlea, angle=-41)
#        od1=dictFlea['od1'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]
#        od1B=dictFlea['od1'][roiB[0]:roiB[1],roiB[2]:roiB[3]]
#        od2=dictFlea['od2'][roiFlea[0]:roiFlea[1],roiFlea[2]:roiFlea[3]]   
#        od2B=dictFlea['od2'][roiB[0]:roiB[1],roiB[2]:roiB[3]]
#        num1=np.sum(od1)-np.sum(od1B)
#        num2=np.sum(od2)-np.sum(od2B)
#        imbal=(num1-num2)/(num1+num2)
#        signalGood=((num1>0) & (num2>0))
#        print num1, num2, imbal, signalGood
#    else:
#        imbal = nan
#        signalGood=False
#    
#    
#    if bgndRidge:
#        odRoi=fringeremoval(filelist,filenum, roi, fileroot=fileroot)
#        
#    else:
#        odRoi=dict1['rotODcorr'][roi[0]:roi[1],roi[2]:roi[3]]   
#    
#    delta=30
#    xM=np.array([])
#    yM=np.array([])
#    fM=np.array([])
#    xvals = np.arange(roi[2]-delta,roi[3]+delta)
#    yvals = np.arange(roi[0]-delta,roi[1]+delta)
#    for x in xvals:
#        for y in yvals:
#            if (x>roi[2] and x<roi[3] and y>roi[0] and y<roi[1])==False:
#                xM=np.append(xM,x)
#                yM=np.append(yM,y)
#                fM=np.append(fM,dict1['rotODcorr'][y,x])
#                
#    (A,B,C), cov = optimize.curve_fit(plane,(xM,yM),fM,p0=(0.03,(roi[3]-roi[2])/2.0,(roi[1]-roi[0]/2.0)))
##    print A,B,C,D,F,G
#    (x,y)=np.meshgrid(np.arange(roi[2],roi[3]),np.arange(roi[0],roi[1]))
#    fitted_plane=plane((x,y),A,B,C).reshape(np.arange(roi[0],roi[1]).size,np.arange(roi[2],roi[3]).size)
#    odRoi1=odRoi-fitted_plane
#    
#    odFiltered=snd.filters.gaussian_filter(odRoi1,1.0)
#    
#    if draw:      
#        fig1=plt.figure()
#        pan1=fig1.add_subplot(1,1,1)
#        pan1.imshow(odFiltered,vmin=-0.15,vmax=0.3)
#        
#    
#    odFiltLin = np.sum(odFiltered,axis=0)
#        
##    peakL=np.max(odFiltLin[:120])
##    xGuessL=np.float( np.where(odFiltLin==peakL)[0])
##    print xGuessL
##    peakR=np.max(odFiltLin[120:])
##
##    xGuessR=np.float( np.where(odFiltLin==peakR)[0])
##    print xGuessR
##    peakT=np.max(odFiltLin)
##    xGuessT=np.float( np.where(odFiltLin==peakT)[0])
##    xGuess2=175
##    xGuess1=xGuess2-55.0
##    xGuess3=xGuess2+52.0
##    
#    if draw:
#        figure2=plt.figure()
#        panel2=figure2.add_subplot(1,1,1)
#        panel2.plot(odFiltLin,'bo')
#    
##    p=(1.5,xGuess1,10.0,peak2,xGuess2,10.0,0.5,xGuess3,10.0,0,0,0)
##    (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3, np.arange(odFiltLin.size), odFiltLin, p0 =p )
##    data_fitted = gaussLin3(np.arange(odFiltLin.size), *(amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C))
#    
##    if draw:
##        panel2.plot(data_fitted,'b-')
#  
##    xCent=(xGuessR+xGuessL)/2.0
#
#    peak=np.max(odFiltLin[160:220])
#    xGuess=np.float(np.where(odFiltLin==peak)[0])
#    xCent=xGuess    
#    
#
#
##    xGuess1=xGuess2-55.0
##    xGuess3=xGuess2+54.0
#
#        
##    p=(1.0,10.0,1.0,xGuess2,10.0,1.0,10.0,4.0,0,0)
##    (amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3const, np.arange(odFiltLin2.size), odFiltLin2,p0=p )
##    data_fitted = gaussLin3const(np.arange(odFiltLin2.size), *(amp1,sigma1,amp2,x02,sigma2,amp3,sigma3,A,B,C))
#
#    
#   # yCent=0.865*xGuess2+54 #odtkick negative
#  #  yCent=1.32*xGuess2-3.3# odtkick positive
#  #  yCent=-0.055*xGuess2**2.0+13.78*xGuess2-690.372
#
#
#    
#        
#    #for rf only:
#    odFiltVert = np.sum(odFiltered,axis=1)  
#    plot(odFiltVert)
#    peakRf=np.max(odFiltVert)
#    yCent = np.float( np.where(odFiltVert==peakRf)[0])
#    print yCent
#    print xCent
##    xCent=165
#    
#  #  print xCent,xGuess2, yCent
##    print yCent
##    yCent=135
#    
#    norm=mpl.colors.Normalize(vmin=-0.15,vmax=0.3)
#    im = Image.fromarray(np.uint8(plt.cm.jet(norm(odFiltered))*255))
#    y1=63
#    x1=58
#    yrf=2
# #   offsets = np.array([[0,0],[0,y1],[0,-y1],[0,y1*2],[0,-y1*2],[-x1,-y2],[-x1,-y2+y1],[-x1,-y2-y1],[-x1,-y2+y1*2],[-x1,-y2-y1*2],[x1,y2],[x1,y2+y1],[x1,y2-y1],[x1,y2+y1*2],[x1,y2-y1*2]])#np.array([[0,0],[0,69],[0,-69],[-58,49],[-58,119],[-58,-21],[57,-47],[57,21],[57,-115]])
# #   offsetsShort=   np.array([[0,0],[0,y1],[0,-y1],[-x1,-y2],[-x1,-y2+y1],[-x1,-y2-y1],[x1,y2],[x1,y2+y1],[x1,y2-y1]]) 
#    rfOffsets = np.array([[0,0],[0,y1],[0,-y1],[-x1,yrf],[-x1,y1+yrf],[-x1,-y1+yrf],[x1,-yrf],[x1,y1-yrf],[x1,-y1-yrf],[-2*x1,yrf],[-2*x1,y1+yrf],[-2*x1,-y1+yrf],[2*x1,-yrf],[2*x1,y1-yrf],[2*x1,-y1-yrf]])    
#    rfOffsetsShortest=np.array([[0,0],[-x1,yrf],[x1,-yrf],[-2*x1,yrf],[2*x1,-yrf]])
#    offsets=rfOffsetsShortest  
#  #  offsets=offsetsShort
#    counts=np.zeros(offsets.shape[0])    
#    for i in np.arange(offsets.shape[0]):
#        offs=offsets[i]
#        counts[i]= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=17,draw=False)[0]
##    print counts
#    maxInd=np.where(counts==np.max(counts))[0][0]
#    
#    allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],r=17)
##    print allcounts
#    i=0
#    while ((np.abs(allcounts[4]-allcounts[3])>2.0) & (i<20)):
#        if (allcounts[4]-allcounts[3])>0:
#            yCent=yCent+1
#       #     print "new yCent = " +str(yCent)
#            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
#        else: 
#            yCent=yCent-1
#      #      print "new yCent = " +str(yCent)
#            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
#        i=i+1
##    print i
#    i=0
#    while ((np.abs(allcounts[2]-allcounts[1])>2.0) & (i<20)):
#        if (allcounts[2]-allcounts[1])>0:
#            xCent=xCent+1
#     #       print "new xCent = " +str(xCent)blochOscFractionsV2(fileroot,filelist,roi,key,plot=True,xlabel='',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea
#            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
#        else: 
#            xCent=xCent-1
#     #       print "new xCent = " +str(xCent)a[6]['Note']
#            allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
#        i=i+1
##    print i
# #   print allcounts
#    
#    for i in np.arange(offsets.shape[0]):
#        offs=offsets[i]
#        r=20.0   
#        if i==0:
#            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
##            if count>0:
##                count -= 6.0
#          #  bgnd1=getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1]+2.0*r,r=r, draw=draw,color=(0,0,0))[0]
#            bgnd=getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1]-2.0*r,r=r, draw=draw,color=(0,0,0))[0]
#            counts[i]=count-bgnd
#        if np.mod(i,offsets.size/3)==1:
#            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
##            if count>0:
##                count -= 6.0
#            bgnd=getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1]-2.0*r,r=r, draw=draw,color=(0,0,0))[0]
#            counts[i]=count-bgnd
#        else:
#            count= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=r, draw=draw)[0]
##            if count>0:
##                count -= 6.0
#            bgnd=getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1]-2.0*r,r=r, draw=draw,color=(0,0,0))[0]
#            counts[i]=count-bgnd
#    if draw:
#        fig4=plt.figure()
#        pan4=fig4.add_subplot(1,1,1)
#        pan4.imshow(im)
#        pan4.set_title(filename)
#        
#    return counts, odFiltered, xCent,yCent, imbal,signalGood,  dict1
#    
#def blochOscFractionsV2(fileroot,filelist,roi,key,plot=True,xlabel='',checkField=True,filerootFlea=filerootFlea,roiFlea=roiFlea):
#   # filerange = np.arange(filestart,filestop+1)    
#    
#    num1=np.zeros(filelist.size)
#    num2=np.zeros(filelist.size)
#    num3=np.zeros(filelist.size)
#    fractionP=np.zeros(filelist.size)
#    fraction0=np.zeros(filelist.size)
#    fractionM=np.zeros(filelist.size)
#    imbalArray=np.ones(filelist.size)
#    qlist=np.zeros(filelist.size)
#    signalGood=np.ones(filelist.size,dtype=bool)
#   # unindexedWave=np.zeros(filelist.size)
#   # filename = fileroot + "_"+ str(filelist[0]).zfill(4) + ".ibw"
#
#    (counts, odFiltered, xCent,yCent,imbal,signalG, dict1)=blochOscOneFileV2(fileroot,filelist[0], roi, draw=False)
#    infoString=dict1["Note"]
#    waveDict=getIndexedWaves(infoString)
#    waveDict=waveDict.fromkeys(waveDict.keys(),[])
#    odFiltAll=np.zeros((filelist.size,odFiltered.shape[0],odFiltered.shape[1]))
#    xCenters=np.zeros(filelist.size)
#    yCenters=np.zeros(filelist.size)
#    for ind, filenum in enumerate(filelist):
#       # filename = fileroot + "_"+ str(filenum).zfill(4) + ".ibw"
#        #print filename
#        (counts, odFiltered, xCent, yCent,imbal, signalG,dict1)=blochOscOneFileV2(fileroot,filenum, roi, draw=False,checkField=checkField,filerootFlea=filerootFlea,roiFlea=roiFlea)
#        print filenum
#        xCenters[ind]=xCent
#        yCenters[ind]=yCent
#        if checkField:
#            imbalArray[ind]=imbal
#            signalGood[ind]=signalG
#        infoString=dict1["Note"]
#        waveDictLocal=getIndexedWaves(infoString)
#        for k in waveDict.iterkeys():
#            
#            waveDict[k]=np.append(waveDict[k],waveDictLocal[k])
#
#  #      unindexedWave[ind]=Generic.StringByKey('nan', key, infoString, '=', '\n')
#        odFiltAll[ind]=odFiltered
#        
##
##        roiNum=counts.size/3
##        num1[ind]=np.sum(counts[roiNum:roiNum*2])
##        if num1[ind]<0.0:
##            print 'num1 < 0 in file ' + str(filenum)
##        num2[ind]=np.sum(counts[0:roiNum])
##        if num2[ind]<0.0:
##            print 'num2 < 0 in file ' + str(filenum)
##        num3[ind]=np.sum(counts[roiNum*2:roiNum*3])
##        if num3[ind]<0.0:
##            print 'num3 < 0 in file ' + str(filenum)
#        num1[ind]=counts[0]
#        num2[ind]=0.0#counts[0]
#        num3[ind]=counts[1]
#        total=num1[ind]+num2[ind]+num3[ind]
#        fractionP[ind]=num1[ind]/total
#        fraction0[ind]=num2[ind]/total
#        fractionM[ind]=num3[ind]/total
#        print num1[ind], num2[ind], num3[ind]
#        print fractionP[ind],fraction0[ind],fractionM[ind]
#        
#    if plot:
#        figure3=plt.figure()
#        panel3=figure3.add_subplot(1,1,1)
#        if key == 'ind':
#            xvals=range(filelist.size)
#        else:
#            xvals =waveDict[key]
#            
#        if xlabel=='':
#            if key == 'ind':
#                xlabel = 'shot number'
#            else:
#                xlabel = key
#        panel3.plot(xvals,fractionP,'bo', label='mF=+1')
#        panel3.plot(xvals,fraction0,'go', label='mF=0')
#        panel3.plot(xvals,fractionM,'ro', label='mF=-1')
#        panel3.set_xlabel(xlabel)
#        panel3.set_ylabel('Fractional population')
#        legend()
#        
#        figure4=plt.figure()
#        pan4=figure4.add_subplot(1,1,1)
#        pan4.plot(xvals, qlist,'bo')
#        pan4.set_ylabel(r'Quasimomentum [$k_L$]')
#        pan4.set_xlabel('Oscillation time [s]')
#        
#        figure = plt.figure()
#        gs = gridspec.GridSpec(3,filelist.size)
#        gs.update(left=0.01, right=0.99, wspace=0.005, hspace=0)
#        if key == 'ind':
#            sort=range(filelist.size)
#            tsorted=range(filelist.size)
#        else:
#            sort=np.argsort(waveDict[key])
#            tsorted = waveDict[key][sort]
#
#        w=30
#        for i in range(filelist.size):
#            j=sort[i]
#            panel1=figure.add_subplot(gs[i])
#            panel1.imshow(odFiltAll[j][:,int(xCenters[j]-49-w):int(xCenters[j]-49+w)],vmin=-0.05, vmax=0.3)
#            panel1.yaxis.set_ticks([])  
#            panel1.xaxis.set_ticks([])
#    
#            panel2=figure.add_subplot(gs[i+filelist.size])
#            panel2.imshow(odFiltAll[j][:,int(xCenters[j]-w):int(xCenters[j]+w)],vmin=-0.05, vmax=0.3)
#            panel2.yaxis.set_ticks([])  
#            panel2.xaxis.set_ticks([])
#            
#            panel3=figure.add_subplot(gs[i+2*filelist.size])
#            panel3.imshow(odFiltAll[j][:,int(xCenters[j]+49-w):int(xCenters[j]+49+w)],vmin=-0.05, vmax=0.3)
#            panel3.yaxis.set_ticks([])  
#            panel3.xaxis.set_ticks([])      
#        
#    return num1,num2,num3,fractionP,fraction0,fractionM,waveDict,imbalArray,signalGood, xCenters, yCenters
#        
#def blochOscOneFile(filename, roi, fitBgnd=True):
#    dict1 =readIgor.processIBW(filename, angle=-37)
#    figure=plt.figure()
#    panel=figure.add_subplot(1,1,1)
#    panel.imshow(dict1['rotOD'],vmin=-0.15, vmax=0.3)
#    xmin,xmax=panel.get_xlim()
#    panel.axhspan(roi[0],roi[1],xmin=(roi[2]-xmin)/(xmax-xmin),xmax=(roi[3]-xmin)/(xmax-xmin), ls='dashed', fill=False)
#    
#    infoString=dict1["Note"]
#    waveDict=readIgor.getIndexedWaves(infoString)
#    print waveDict
#    
#    odRoi=dict1['rotOD'][roi[0]:roi[1],roi[2]:roi[3]]
#    
#    if fitBgnd:
#        delta=30
#        xM=np.array([])
#        yM=np.array([])
#        fM=np.array([])
#        xvals = np.arange(roi[2]-delta,roi[3]+delta)
#        yvals = np.arange(roi[0]-delta,roi[1]+delta)
#        for x in xvals:
#            for y in yvals:
#                if (x>roi[2] and x<roi[3] and y>roi[0] and y<roi[1])==False:
#                    xM=np.append(xM,x)
#                    yM=np.append(yM,y)
#                    fM=np.append(fM,dict1['rotOD'][y,x])
#                    
#        (A,B,C), cov = optimize.curve_fit(plane,(xM,yM),fM)
#        
#        (x,y)=np.meshgrid(np.arange(roi[2],roi[3]),np.arange(roi[0],roi[1]))
#        fitted_plane=plane((x,y),A,B,C).reshape(np.arange(roi[0],roi[1]).size,np.arange(roi[2],roi[3]).size)
#        odRoi=odRoi-fitted_plane
#
#
#
#    odFiltered=snd.filters.gaussian_filter(odRoi,1.3)
#    
#    figure2=plt.figure()
#    panel2=figure2.add_subplot(1,1,1)
#    panel2.imshow(odFiltered,vmin=-0.15, vmax=0.3)
#    odFiltLin = np.sum(odFiltered,axis=0)
#    
#    peak2=np.max(odFiltLin[110:170])
#    print peak2
#    xGuess2=np.float( np.where(odFiltLin==peak2)[0])
#    print xGuess2
#
#    xGuess1=xGuess2-55.0
#    xGuess3=xGuess2+59.0
#        
#    figure3=plt.figure()
#    panel3=figure3.add_subplot(1,1,1)
#    panel3.plot(odFiltLin,'bo')
#    
#    p=(1.5,xGuess1,10.0,peak2,xGuess2,10.0,0.5,xGuess3,10.0,0,0,0)
#    (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3, np.arange(odFiltLin.size), odFiltLin, p0 =p )
#    data_fitted = gaussLin3(np.arange(odFiltLin.size), *(amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C))
#    
#    print (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C)
#    print np.diag(covariance)
##    if np.diag(covariance)[2]>1:
##        p=(1.5,28.,4.0,2.5,111.0,15.0,1.0,198.0,10.0,2,-0.005)
##        (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B), covariance = optimize.curve_fit(gaussLin3, np.arange(odFiltLin.size), odFiltLin, p0 =p )
##        data_fitted = gaussLin3(np.arange(odFiltLin.size), *(amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B))
##        print np.diag(covariance)
#
#    panel3.plot(data_fitted,'b-')
#    
#    num1=np.sqrt(np.pi)*amp1*np.abs(sigma1)
#
#    num2=np.sqrt(np.pi)*amp2*np.abs(sigma2)
#
#    num3=np.sqrt(np.pi)*amp3*np.abs(sigma3)
#
#    
#    return num1, num2, num3
#    
#    
#def blochOscFractions(fileroot,filelist,roi,key,plot=True,xlabel='',fitBgnd=True):
#   # filerange = np.arange(filestart,filestop+1)    
#    
#    num1=np.zeros(filelist.size)
#    num2=np.zeros(filelist.size)
#    num3=np.zeros(filelist.size)
#    fraction1=np.zeros(filelist.size)
#    fraction2=np.zeros(filelist.size)
#    fraction3=np.zeros(filelist.size)
#    
#    filename = fileroot + "_"+ str(filelist[0]).zfill(4) + ".ibw"
#    dict1 =readIgor.processIBW(filename,angle=-37)
#    odRoi=dict1['rotOD'][roi[0]:roi[1],roi[2]:roi[3]]
#    infoString=dict1["Note"]
#    waveDict=readIgor.getIndexedWaves(infoString)
#    waveDict=waveDict.fromkeys(waveDict.keys(),[])
#    odFiltered=np.zeros((filelist.size,odRoi.shape[0],odRoi.shape[1]))
#    
#    for ind, filenum in enumerate(filelist):
#        filename = fileroot + "_"+ str(filenum).zfill(4) + ".ibw"
#        dict1 =readIgor.processIBW(filename, angle=-37)
#        
#        infoString=dict1["Note"]
#        waveDictLocal=readIgor.getIndexedWaves(infoString)
#        for k in waveDict.iterkeys():
#            
#            waveDict[k]=np.append(waveDict[k],waveDictLocal[k])
# 
#  
#     
#        odRoi=dict1['rotOD'][roi[0]:roi[1],roi[2]:roi[3]]
#        
#        if fitBgnd:
#            delta=30
#            xM=np.array([])
#            yM=np.array([])
#            fM=np.array([])
#            xvals = np.arange(roi[2]-delta,roi[3]+delta)
#            yvals = np.arange(roi[0]-delta,roi[1]+delta)
#            for x in xvals:
#                for y in yvals:
#                    if (x>roi[2] and x<roi[3] and y>roi[0] and y<roi[1])==False:
#                        xM=np.append(xM,x)
#                        yM=np.append(yM,y)
#                        fM=np.append(fM,dict1['rotOD'][y,x])
#            (A,B,C), cov = optimize.curve_fit(plane,(xM,yM),fM)
#            
#            (x,y)=np.meshgrid(np.arange(roi[2],roi[3]),np.arange(roi[0],roi[1]))
#            fitted_plane=plane((x,y),A,B,C).reshape(np.arange(roi[0],roi[1]).size,np.arange(roi[2],roi[3]).size)
#            odRoi=odRoi-fitted_plane
#
#        odFiltered[ind,:,:]=snd.filters.gaussian_filter(odRoi,1.3)
#        
#        bound1=110
#        bound2=170
#        odFiltLin = np.sum(odFiltered[ind,:,:],axis=0)
#        peak2=np.max(odFiltLin[bound1:bound2])
#        xGuess2=np.float( np.where(odFiltLin==peak2)[0])
#        xGuess1=xGuess2-55.0
#        xGuess3=xGuess2+59.0
#            
#        print filenum
#        p=(1.5,xGuess1,10.0,peak2,xGuess2,10.0,0.5,xGuess3,10.0,0,0,0)
#        (amp1,x01,sigma1,amp2,x02,sigma2,amp3,x03,sigma3,A,B,C), covariance = optimize.curve_fit(gaussLin3, np.arange(odFiltLin.size), odFiltLin, p0 = p)
#
#        
#        num1[ind]=np.sqrt(np.pi)*amp1*np.abs(sigma1)
#        if num1[ind]<0.0:
#            print 'num1 < 0 in file ' + str(filenum)
#        num2[ind]=np.sqrt(np.pi)*amp2*np.abs(sigma2)
#        if num2[ind]<0.0:
#            print 'num2 < 0 in file ' + str(filenum)
#        num3[ind]=np.sqrt(np.pi)*amp3*np.abs(sigma3)
#        if num3[ind]<0.0:
#            print 'num3 < 0 in file ' + str(filenum)
#        total=num1[ind]+num2[ind]+num3[ind]
#        fraction1[ind]=num1[ind]/total
#        fraction2[ind]=num2[ind]/total
#        fraction3[ind]=num3[ind]/total
#        
#    if plot:
#        figure3=plt.figure()
#        panel3=figure3.add_subplot(1,1,1)
#        if key == 'ind':
#            xvals=range(filelist.size)
#        else:
#            xvals =waveDict[key]
#            
#        if xlabel=='':
#            if key == 'ind':
#                xlabel = 'shot number'
#            else:
#                xlabel = key
#        panel3.plot(xvals,fraction1,'bo', label='mF=+1')
#        panel3.plot(xvals,fraction2,'go', label='mF=0')
#        panel3.plot(xvals,fraction3,'ro', label='mF=-1')
#        panel3.set_xlabel(xlabel)
#        panel3.set_ylabel('Fractional population')
#        legend()
#        
#        figure = plt.figure()
#        gs = gridspec.GridSpec(3,filelist.size)
#        if key == 'ind':
#            sort=range(filelist.size)
#            tsorted=range(filelist.size)
#        else:
#            sort=np.argsort(waveDict[key])
#            tsorted = waveDict[key][sort]
#
#        
#        for i in range(filelist.size):
#            j=sort[i]
#            panel1=figure.add_subplot(gs[i])
#            panel1.imshow(odFiltered[j,:,:bound1],vmin=-0.05, vmax=0.3)
#            panel1.yaxis.set_ticks([])  
#            panel1.xaxis.set_ticks([])
#    
#            panel2=figure.add_subplot(gs[i+filelist.size])
#            panel2.imshow(odFiltered[j,:,bound1:bound2],vmin=-0.05, vmax=0.3)
#            panel2.yaxis.set_ticks([])  
#            panel2.xaxis.set_ticks([])
#            
#            panel3=figure.add_subplot(gs[i+2*filelist.size])
#            panel3.imshow(odFiltered[j,:,bound2:],vmin=-0.05, vmax=0.3)
#            panel3.yaxis.set_ticks([])  
#            panel3.xaxis.set_ticks([])        
#    return num1,num2,num3,fraction1,fraction2,fraction3,waveDict
