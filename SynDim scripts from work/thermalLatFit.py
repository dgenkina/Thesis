# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:48:04 2017

@author: dng5
"""

import numpy as np
import scipy.ndimage as snd
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import IgorBin
import time

fileroot = 'Y:/Data/2017/October/18/PIXIS_18Oct2017' 
angle=-39

PIXIS_background_filename ='Y:/Data/2017/October/13/PIXIS_13Oct2017_10003.ibw'
bgndImage=IgorBin.LoadIBW(PIXIS_background_filename)['Data']
b1=np.array(bgndImage[:,0:bgndImage.shape[1]/2], dtype=float)
b2=np.array(bgndImage[:,bgndImage.shape[1]/2:], dtype=float)

filelist=np.arange(4,34)

roi=np.array([430, 690, 300,650])
key = 'pulseDelay'
fitProbes=True
filelistProbes=np.arange(4,34)

def getProbeReconstruction(fileroot, filelistProbes, roi):
    filename=fileroot+"_"+ str(filelistProbes[0]).zfill(4) + ".ibw"
    dict1 =processIBW(filename, angle=angle,bgnd1=b1,bgnd2=b2)  
    rotProbe=snd.interpolation.rotate(dict1['Raw2'],angle)[roi[0]:roi[1],roi[2]:roi[3]]
    N=rotProbe.flatten().size
    print 'number of pixels in roi = ' + str(N)
    t1=time.clock()
    probeMatrix=np.ones((N,filelistProbes.size+1))
    for ind, filenum in enumerate(filelistProbes):
        filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
        dict1 =processIBW(filename, angle=angle,bgnd1=b1,bgnd2=b2)  
        rotProbe=dict1['rotRaw2'][roi[0]:roi[1],roi[2]:roi[3]]
        probeMatrix[:,ind+1]=rotProbe.flatten()
    probeMatrixT=probeMatrix.transpose()
    t2=time.clock()
    print 'constructed probe matrix in '+str(t2-t1)
    return probeMatrix, probeMatrixT

'''Calculate matrix of probes'''
if fitProbes:
    probeMatrix,probeMatrixT=getProbeReconstruction(fileroot, filelistProbes, roi)
else:
    probeMatrix=1.0
    probeMatrixT=1.0


def thermalOneFile(filenum, plot=True):
    filename = fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
    print filename
    dict1=processIBW(filename, angle=angle,bgnd1=b1,bgnd2=b2)

    
    odRoi=dict1['rotODcorr'][roi[0]:roi[1],roi[2]:roi[3]]   
    atomsRoi=dict1['rotRaw1'][roi[0]:roi[1],roi[2]:roi[3]]
        
    if fitProbes:
        weightArray=np.ones(odRoi.shape)
        weightArray[70:250,120:240] = 0


        xTw=probeMatrixT*weightArray.flatten()
        rhs=np.dot(xTw,atomsRoi.flatten())
        lhs=np.dot(xTw,probeMatrix)
        beta=np.linalg.solve(lhs,rhs)
        newProbe=np.dot(probeMatrix,beta).reshape(atomsRoi.shape)
        newOd=-np.log(atomsRoi/newProbe)
        newOd = zeroNansInfsVector(newOd)
        odRoi = newOd+(newProbe-atomsRoi)/(IsatCounts*imagingTime)


    
    odLin = np.sum(odRoi[:,160:200],axis=1)
    odLinNorm=odLin/np.sum(odLin)
    if plot:
        figure=plt.figure()
        pan=figure.add_subplot(111)
        pan.imshow(odRoi)
        
        if fitProbes:
            figure=plt.figure()
            pan=figure.add_subplot(111)
            pan.imshow(newOd)
        
        fig2=plt.figure()
        pan2=fig2.add_subplot(111)
        pan2.plot(odLin,'bo')
    
    return dict1, odLin, odLinNorm
    
    
dict1, odLin, odLinNorm = thermalOneFile(filelist[0], plot=False)   
infoString=dict1["Note"]
waveDict=getIndexedWaves(infoString)
waveDict=waveDict.fromkeys(waveDict.keys(),[])    
buildOD=np.zeros((roi[1]-roi[0],filelist.size))

for i,filenum in enumerate(filelist):
    dict1, odLin, odLinNorm = thermalOneFile(filenum, plot=False)
   
    buildOD[:,i] = odLinNorm
    
    infoString=dict1["Note"]
    waveDictLocal=getIndexedWaves(infoString)
    for k in waveDict.iterkeys():
        
        waveDict[k]=np.append(waveDict[k],waveDictLocal[k])
        
sort = np.argsort(waveDict[key])
tlist=waveDict[key][sort]
figure=plt.figure()
pan=figure.add_subplot(111)
pan.imshow(buildOD[:,sort])

slicewidth=2
test=buildOD.reshape(buildOD.shape[0]/slicewidth,slicewidth,buildOD.shape[1])
test2=np.sum(test,axis=1)

figure=plt.figure()
pan=figure.add_subplot(111)
pan.imshow(test2[:,sort], aspect=0.000004,extent=(waveDict[key][sort][0],waveDict[key][sort][-1],roi[1]/slicewidth,roi[2]/slicewidth))
   
flist=np.zeros(buildOD.shape[0])
uncertFlist=  np.zeros(buildOD.shape[0]) 
fitGood=np.ones(buildOD.shape[0], dtype=bool)
for i in np.arange(buildOD.shape[0]):
    try:
        A,flist[i],phi,offset,resids,chi2,uncerts = sineFit(tlist,buildOD[i,sort],'lattice pulse time','population',p0=(-0.001,500,3.9,0),plot=False)
        print chi2,uncerts
        uncertFlist[i]=uncerts[1]
        if uncerts[1]>80.0:
            fitGood[i]=False
        elif np.isnan(uncerts[1]):
            fitGood[i]=False
    except RuntimeError:
        print 'Unable to fit line ' + str(i)
        fitGood[i]=False

easyPlot(np.arange(buildOD.shape[0])[fitGood],np.abs(flist[fitGood]))

    