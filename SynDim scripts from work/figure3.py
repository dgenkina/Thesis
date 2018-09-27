# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:21:44 2017

@author: dng5
"""


import matplotlib.pyplot as plt
#import numpy as np



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np
import lineFit2

# Use me to produce pdf files
# from matplotlib.backends.backend_pgf import FigureCanvasPgf
# matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8


rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['text.latex.preamble'] = [r'\usepackage{cmbright}',r'\usepackage{amsmath}']


rcParams['axes.linewidth'] = 0.75
rcParams['lines.linewidth'] = 0.75

rcParams['xtick.major.size'] = 3      # major tick size in points
rcParams['xtick.minor.size'] = 2      # minor tick size in points
rcParams['xtick.major.width'] = 0.75       # major tick width in points
rcParams['xtick.minor.width'] = 0.75      # minor tick width in points

rcParams['ytick.major.size'] = 3      # major tick size in points
rcParams['ytick.minor.size'] = 2      # minor tick size in points
rcParams['ytick.major.width'] = 0.75       # major tick width in points
rcParams['ytick.minor.width'] = 0.75      # minor tick width in points

ColorMap = "afmhot"

datafileF1p=np.load('08Mar2017_files_154-183.npz')
datafileF1m=np.load('08Mar2017_files_64-93.npz')
datafileF2m=np.load('28Feb2017_files_244-273.npz')
datafileF2p=np.load('28Feb2017_files_214-243.npz')

cdictC = {'red':   [(0.0,  1.0, 1.0), 
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.5, 0.5)],

         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  0.5, 0.5)]}
                   
cyans = mpl.colors.LinearSegmentedColormap('Cyans', cdictC)
plt.register_cmap(cmap=cyans)        

cdictM = {'red':   [(0.0,  1.0, 1.0), 
                   (1.0,  0.5, 0.5)],

         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,   0.5,  0.5)]}
                   
magentas = mpl.colors.LinearSegmentedColormap('Magentas', cdictM)
plt.register_cmap(cmap=magentas)

#pickFrac1p=(np.mod(np.arange(datafileF1p['xCenters'].size),2)==0)
#odF1p=datafileF1p['odFiltAll'][pickFrac1p]
#xCentersF1p=datafileF1p['xCenters'][pickFrac1p]
#tlistF1p=datafileF1p['tlist'][pickFrac1p]
#pickFrac1m=(np.mod(np.arange(datafileF1m['xCenters'].size),2)==0)
#odF1m=datafileF1m['odFiltAll'][pickFrac1m]
#xCentersF1m=datafileF1m['xCenters'][pickFrac1m]
#tlistF1m=datafileF1m['tlist'][pickFrac1m]

def demap(odIn,y0,kL,bC,klist):
    odOut=np.zeros((2*kL,odIn.shape[1]))
    
    for i in range(odOut.shape[0]):
        for j in range(odOut.shape[1]):
            norm=0
            for m in np.array([-4,-2,0,2,4]):
                bCind=np.where(np.abs(klist-(i-kL+m*kL)/kL)==np.min(np.abs(klist-(i-kL+m*kL)/kL)))
                try:
                    odOut[i,j]=odOut[i,j]+odIn[i-kL+y0+m*kL,j]#/bC[bCind]
                except IndexError:
                    odOut[i,j]+=0.0
                norm=1.0#+=1.0/(bC[bCind]**2.0)
            odOut[i,j]=odOut[i,j]/norm
    return odOut

pickFrac2p=(np.mod(np.arange(datafileF2p['xCenters'].size),1)==0)
odF2p=datafileF2p['odFiltAll'][pickFrac2p]
xCentersF2p=datafileF2p['xCenters'][pickFrac2p]
tlistF2p=datafileF2p['tlist'][pickFrac2p]
qlistF2p=datafileF2p['qlist'][pickFrac2p]
pickFrac2m=(np.mod(np.arange(datafileF2m['xCenters'].size),1)==0)
odF2m=datafileF2m['odFiltAll'][pickFrac2m]
xCentersF2m=datafileF2m['xCenters'][pickFrac2m]
tlistF2m=datafileF2m['tlist'][pickFrac2m]
qlistF2m=datafileF2m['qlist'][pickFrac2m]



#sort1p=np.argsort(tlistF1p)
#tsorted1p = tlistF1p[sort1p]
#w=30
#for i in range(tsorted1p.size):
#    j=sort1p[i]
#
#    panel2=figure.add_subplot(gs[i+xCentersF1m.size+offset])
#    panel2.imshow(odF1p[j][:,int(xCentersF1p[j]-49-w):int(xCentersF1p[j]-49+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel2.yaxis.set_ticks([])  
#    panel2.xaxis.set_ticks([])
#    panel2.set_frame_on(False)
#
#    panel3=figure.add_subplot(gs[i+xCentersF1m.size+numTot1+offset])
#    panel3.imshow(odF1p[j][:,int(xCentersF1p[j]-w):int(xCentersF1p[j]+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel3.yaxis.set_ticks([])  
#    panel3.xaxis.set_ticks([])
#    panel3.set_frame_on(False)
#    
#    panel4=figure.add_subplot(gs[i+xCentersF1m.size+2*numTot1+offset])
#    panel4.imshow(odF1p[j][:,int(xCentersF1p[j]+49-w):int(xCentersF1p[j]+49+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel4.yaxis.set_ticks([])  
#    panel4.xaxis.set_ticks([])  
#    panel4.set_frame_on(False)
#
#sort1m=np.argsort(tlistF1m)
#tsorted1m = tlistF1m[sort1m]
w=30
#for i in range(tsorted1m.size):
#    j=sort1m[i]
#
#    panel2=figure.add_subplot(gs[xCentersF1p.size-i+offset])
#    panel2.imshow(odF1m[j][:,int(xCentersF1m[j]-49-w):int(xCentersF1m[j]-49+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel2.yaxis.set_ticks([])  
#    panel2.xaxis.set_ticks([])
#    panel2.set_frame_on(FdatafileF2p['tlist']alse)
##    if i==tsorted1m.size-1:
##        panel2.set_ylabel(r'$m_F= +1$',labelpad=22,rotation='horizontal',va='center')
#
#    panel3=figure.add_subplot(gs[xCentersF1p.size+numTot1-i+offset])
#    panel3.imshow(odF1m[j][:,int(xCentersF1m[j]-w):int(xCentersF1m[j]+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel3.yaxis.set_ticks([])  
#    panel3.xaxis.set_ticks([])
#    panel3.set_frame_on(False)
##    if i==tsorted1m.size-1:
##        panel3.set_ylabel(r'$m_F= 0$',labelpad=22,rotation='horizontal',va='center')
#        
#    panel4=figure.add_subplot(gs[xCentersF1p.size+2*numTot1-i+offset])
#    panel4.imshow(odF1m[j][:,int(xCentersF1m[j]+49-w):int(xCentersF1m[j]+49+w)],vmin=-0.05, vmax=0.3,interpolation='none')
#    panel4.yaxis.set_ticks([])  
#    panel4.xaxis.set_ticks([])
#    panel4.set_frame_on(False)
##    if i==tsorted1m.size-1:
##        panel4.set_ylabel(r'$m_F= -1$',labelpad=22,rotation='horizontal',va='center')



w=40  
dist=49  
kL=31
y0=175  
file1=np.load('blochCollect_U_4.40.npz')
klist=file1['klist']
bC=file1['bC']

odF2=np.zeros((tlistF2p.size+tlistF2m.size,5,odF2p.shape[1],w*2))
odF2demap=np.zeros((tlistF2p.size+tlistF2m.size,5,kL*2,w*2))

sort2p=np.argsort(tlistF2p)
tsorted2p = tlistF2p[sort2p]  
qsorted2p = qlistF2p[sort2p]  
for i,q in enumerate(qlistF2p):
     if q>-0.8:
         qlistF2p[i] = qlistF2p[i] - 2.0
A,B,dA,dB = lineFit2.lineFit(tlistF2p,qlistF2p,'','')
qlistF2pFit = A*tlistF2p[sort2p]

sort2m=np.argsort(tlistF2m)
tsorted2m = tlistF2m[sort2m] 
for i,q in enumerate(qlistF2m):
     if q<-1.5:
         qlistF2m[i] = qlistF2m[i] + 2.0

A,B,dA,dB = lineFit2.lineFit(tlistF2m,qlistF2m,'','')
 
qlistF2mFit = A*tlistF2m[sort2m]
 
qsorted2m = qlistF2m[sort2m]    

figure = plt.figure()
figure.clear()
figure.set_size_inches(3.3,1.1)

fracnum=6
indlistp=np.arange(0,tsorted2p.size,fracnum)
indlistm=np.arange(0,tsorted2m.size,fracnum)
numTot=indlistp.size+indlistm.size
offset=0
numTot=numTot+offset
gs = gridspec.GridSpec(5,numTot)

gs.update(left=0.0, right=1.0,top=1.0,bottom=0.0, wspace=0.1, hspace=0.0)


vmax=0.4
for i in indlistp:
    j=sort2p[i]
    i=int(i/fracnum)
    panel1=figure.add_subplot(gs[i+indlistm.size+offset])#+4*numTot1
    odF2[j][0]=odF2p[j][:,int(xCentersF2p[j]-2*dist-w):int(xCentersF2p[j]-2*dist+w)]
 #   odF2demap[j][0]=demap(odF2[j][0],y0,kL,bC,klist)
    panel1.imshow(odF2[j][0][y0-kL:y0+kL,15:60].transpose(),cmap='Cyans',vmin=-0.05, vmax=vmax,interpolation='none')
    panel1.yaxis.set_ticks([])  
    panel1.xaxis.set_ticks([])        
    panel1.set_frame_on(False)    
    
    panel2=figure.add_subplot(gs[i+indlistm.size+numTot+offset])
    odF2[j][1]=odF2p[j][:,int(xCentersF2p[j]-dist-w):int(xCentersF2p[j]-dist+w)]
  #  odF2demap[j][1]=demap(odF2[j][1],y0,kL,bC,klist)
    panel2.imshow(odF2[j][1][y0-kL:y0+kL,15:60].transpose(),cmap='Blues',vmin=-0.05, vmax=vmax,interpolation='none')
    panel2.yaxis.set_ticks([])  
    panel2.xaxis.set_ticks([])
    panel2.set_frame_on(False)

    panel3=figure.add_subplot(gs[i+indlistm.size+2*numTot+offset])
    odF2[j][2]=odF2p[j][:,int(xCentersF2p[j]-w):int(xCentersF2p[j]+w)]
 #   odF2demap[j][2]=demap(odF2[j][2],y0,kL,bC,klist)
    panel3.imshow(odF2[j][2][y0-kL:y0+kL,15:60].transpose(),cmap='Greens',vmin=-0.05, vmax=vmax,interpolation='none')
    panel3.yaxis.set_ticks([])  
    panel3.xaxis.set_ticks([])
    panel3.set_frame_on(False)
    
    panel4=figure.add_subplot(gs[i+indlistm.size+3*numTot+offset])
    odF2[j][3]=odF2p[j][:,int(xCentersF2p[j]-w+dist):int(xCentersF2p[j]+w+dist)]
 #   odF2demap[j][3]=demap(odF2[j][3],y0,kL,bC,klist)
    panel4.imshow(odF2[j][3][y0-kL:y0+kL,15:60].transpose(),cmap='Reds',vmin=-0.05, vmax=vmax,interpolation='none')    
    panel4.yaxis.set_ticks([])  
    panel4.xaxis.set_ticks([])  
    panel4.set_frame_on(False)
    
    panel5=figure.add_subplot(gs[i+indlistm.size+4*numTot+offset])
    odF2[j][4]=odF2p[j][:,int(xCentersF2p[j]-w+2*dist):int(xCentersF2p[j]+w+2*dist)]
 #   odF2demap[j][4]=demap(odF2[j][4],y0,kL,bC,klist)
    panel5.imshow(odF2[j][4][y0-kL:y0+kL,15:60].transpose(),cmap='Magentas',vmin=-0.05, vmax=vmax,interpolation='none')
    panel5.yaxis.set_ticks([])  
    panel5.xaxis.set_ticks([]) 
    panel5.set_frame_on(False)


for i in indlistm:
    j=sort2m[i]
    i=int(i/fracnum)
    panel1=figure.add_subplot(gs[indlistm.size-i+offset])#+4*numTot1
    odF2[j+tlistF2p.size][0]=odF2m[j][:,int(xCentersF2m[j]-2*dist-w):int(xCentersF2m[j]-2*dist+w)]
#    odF2demap[j+tlistF2p.size][0]=demap(odF2[j+tlistF2p.size][0],y0,kL,bC,klist)
    panel1.imshow(odF2[j+tlistF2p.size][0][y0-kL:y0+kL,15:60].transpose(),cmap='Cyans',vmin=-0.05, vmax=vmax,interpolation='none')
    panel1.yaxis.set_ticks([])  
    panel1.xaxis.set_ticks([])    
    panel1.set_frame_on(False)
#    if i==tsorted1m.size-1:
#        panel1.set_ylabel(r'$m_F= -2$',labelpad=22,rotation='horizontal',va='center')        
    
    panel2=figure.add_subplot(gs[indlistm.size+numTot-i+offset])
    odF2[j+tlistF2p.size][1]=odF2m[j][:,int(xCentersF2m[j]-dist-w):int(xCentersF2m[j]-dist+w)]
 #   odF2demap[j+tlistF2p.size][1]=demap(odF2[j+tlistF2p.size][1],y0,kL,bC,klist)
    panel2.imshow(odF2[j+tlistF2p.size][1][y0-kL:y0+kL,15:60].transpose(),cmap='Blues',vmin=-0.05, vmax=vmax,interpolation='none')
    panel2.yaxis.set_ticks([])  
    panel2.xaxis.set_ticks([])
    panel2.set_frame_on(False)
#    if i==tsorted1m.size-1:
#        panel2.set_ylabel(r'$m_F= -1$',labelpad=22,rotation='horizontal',va='center')

    panel3=figure.add_subplot(gs[indlistm.size+2*numTot-i+offset])
    odF2[j+tlistF2p.size][2]=odF2m[j][:,int(xCentersF2m[j]-w):int(xCentersF2m[j]+w)]
#    odF2demap[j+tlistF2p.size][2]=demap(odF2[j+tlistF2p.size][2],y0,kL,bC,klist)
    panel3.imshow(odF2[j+tlistF2p.size][2][y0-kL:y0+kL,15:60].transpose(),cmap='Greens',vmin=-0.05, vmax=vmax,interpolation='none')
    panel3.yaxis.set_ticks([])  
    panel3.xaxis.set_ticks([])
    panel3.set_frame_on(False)
#    if i==tsorted1m.size-1:
#        panel3.set_ylabel(r'$m_F= 0$',labelpad=22,rotation='horizontal',va='center')
    
    panel4=figure.add_subplot(gs[indlistm.size+3*numTot-i+offset])
    odF2[j+tlistF2p.size][3]=odF2m[j][:,int(xCentersF2m[j]+dist-w):int(xCentersF2m[j]+dist+w)]
#    odF2demap[j+tlistF2p.size][3]=demap(odF2[j+tlistF2p.size][3],y0,kL,bC,klist)
    panel4.imshow(odF2[j+tlistF2p.size][3][y0-kL:y0+kL,15:60].transpose(),cmap='Reds',vmin=-0.05, vmax=vmax,interpolation='none')
    panel4.yaxis.set_ticks([])  
    panel4.xaxis.set_ticks([]) 
    panel4.set_frame_on(False)
#    if i==tsorted1m.size-1:
#        panel4.set_ylabel(r'$m_F= +1$',labelpad=22,rotation='horizontal',va='center')
    
    panel5=figure.add_subplot(gs[indlistm.size+4*numTot-i+offset])
    odF2[j+tlistF2p.size][4]=odF2m[j][:,int(xCentersF2m[j]+2*dist-w):int(xCentersF2m[j]+2*dist+w)]
#    odF2demap[j+tlistF2p.size][4]=demap(odF2[j+tlistF2p.size][4],y0,kL,bC,klist)
    panel5.imshow(odF2[j+tlistF2p.size][4][y0-kL:y0+kL,15:60].transpose(),cmap='Magentas',vmin=-0.05, vmax=vmax,interpolation='none')
    panel5.yaxis.set_ticks([])  
    panel5.xaxis.set_ticks([]) 
    panel5.set_frame_on(False)
#    if i==tsorted1m.size-1:
#        panel5.set_ylabel(r'$m_F= +2$',labelpad=22,rotation='horizontal',va='center')
    
figure.show()
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure2v7c.jpg', dpi=500)

np.savez('dataForImagesFigure',odF2=np.append(odF2m[sort2m][indlistm][::-1],odF2p[sort2p][indlistp][1:],axis=0),q=np.append(-qlistF2mFit[indlistm][::-1],-qlistF2pFit[indlistp][1:]))