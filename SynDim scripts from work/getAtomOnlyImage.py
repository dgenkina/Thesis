# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:32:06 2016

@author: dng5
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import numpy as np
import scipy.ndimage as snd
from scipy import optimize
import matplotlib.gridspec as gridspec

def plane((x,y), A,B,C):
    pl=A*x+B*y+C
    return pl.ravel()

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


roi=np.array([380, 690, 350,700])
fileroot = 'X:/2016/November/10/PIXIS_10Nov2016'  
filenum=173



filename=fileroot+"_"+ str(filenum).zfill(4) + ".ibw"
dict1 =processIBW(filename, angle=-43.5)

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

odFiltLin = np.sum(odFiltered,axis=0)

peak=np.max(odFiltLin[160:220])
xGuess=np.float(np.where(odFiltLin==peak)[0])
xCent=xGuess 


odFiltVert = np.sum(odFiltered,axis=1)  
peakRf=np.max(odFiltVert)
yCent = np.float( np.where(odFiltVert==peakRf)[0])

norm=mpl.colors.Normalize(vmin=-0.15,vmax=0.3)
im = Image.fromarray(np.uint8(plt.cm.jet(norm(odFiltered))*255))

y1=71
x1=58
yrf=0
rfOffsets = np.array([[0,0],[0,y1],[0,-y1],[-x1,yrf],[-x1,y1+yrf],[-x1,-y1+yrf],[x1,-yrf],[x1,y1-yrf],[x1,-y1-yrf],[-2*x1,yrf],[-2*x1,y1+yrf],[-2*x1,-y1+yrf],[2*x1,-yrf],[2*x1,y1-yrf],[2*x1,-y1-yrf]])    
offsets=rfOffsets

counts=np.zeros(offsets.shape[0])    
for i in np.arange(offsets.shape[0]):
    offs=offsets[i]
    counts[i]= getRoi(odFiltered, im, xCent+offs[0],yCent+offs[1],r=17,draw=False)[0]

maxInd=np.where(counts==np.max(counts))[0][0]

allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],r=17)

i=0
while ((np.abs(allcounts[4]-allcounts[3])>2.0) & (i<20)):
    if (allcounts[4]-allcounts[3])>0:
        yCent=yCent+1
        allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
    else: 
        yCent=yCent-1
        allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
    i=i+1
i=0
while ((np.abs(allcounts[2]-allcounts[1])>2.0) & (i<20)):
    if (allcounts[2]-allcounts[1])>0:
        xCent=xCent+1
        allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
    else: 
        xCent=xCent-1
        allcounts=getRoi(odFiltered, im, xCent+offsets[maxInd,0],yCent+offsets[maxInd,1],draw=False)
    i=i+1

eps=5.0
r=10.0

buildOD=np.zeros((offsets.shape[0],odFiltered.shape[0],odFiltered.shape[1]))
ylen=odFiltered.shape[0]
xlen=odFiltered.shape[1]
y,x = np.ogrid[:ylen, :xlen]
for i in range(offsets.shape[0]):

    mask = np.sqrt((x-xCent-offsets[i,0]-eps)**2.0 + (y-yCent-offsets[i,1])**2.0) +np.sqrt((x-xCent-offsets[i,0]+eps)**2.0 + (y-yCent-offsets[i,1])**2.0) <= 2.0*r
    buildOD[i]=odFiltered
    buildOD[i][np.invert(mask)]=0
    
buildOD=np.sum(buildOD,axis=0)

figure=plt.figure()
pan=figure.add_subplot(1,1,1)
pan.imshow(buildOD)    
    