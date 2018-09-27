# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 16:55:01 2018

@author: dng5
"""

import matplotlib.pyplot as plt
import lineFit2



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np



rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8

rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']


rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}',r'\usepackage{amsmath}']


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
datafileF2p=np.load('28Feb2017_files_244-273.npz')

datafileF2m=np.load('28Feb2017_files_454-483.npz')


#odF2p=datafileF2p['odFiltAll']
#xCentersF2p=datafileF2p['xCenters']
#tlistF2p=datafileF2p['tlist']
#qlistF2p=datafileF2p['qlist']
#for i,q in enumerate(qlistF2p):
#    if q<-1.5:
#        qlistF2p[i] = qlistF2p[i] + 2.0
#        
#A,B,dA,dB = lineFit2.lineFit(tlistF2p,qlistF2p,'','')
#
#qlistF2pFit = tlistF2p*A
#sortF2p=np.argsort(qlistF2pFit)


odF2m=datafileF2m['odFiltAll']
xCentersF2m=datafileF2m['xCenters']
tlistF2m=datafileF2m['tlist']
qlistF2m=datafileF2m['qlist']
for i,q in enumerate(qlistF2m):
    if q>-0.8:
        qlistF2m[i] = qlistF2m[i] - 2.0
        
A,B,dA,dB = lineFit2.lineFit(tlistF2m,qlistF2m,'','')

qlistF2mFit = tlistF2m*A
sortF2m=np.argsort(qlistF2mFit)

        
figure = plt.figure()
figure.clear()
figure.set_size_inches(3.0,1.7)

gs = gridspec.GridSpec(5,1)

gs.update(left=0.05, right=0.95, top=1.0, bottom = 0.2,wspace=0,hspace=0)

vmax=0.4
vmin=0.0

pan0=figure.add_subplot(gs[0])
pan0.imshow(odF2m[sortF2m][-1].transpose()[45:92,85:270],cmap='Cyans',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
pan0.yaxis.set_ticks([14]) 
pan0.set_yticklabels([]) 
pan0.xaxis.set_ticks([])        
pan0.spines["bottom"].set_visible(False) 
pan1=figure.add_subplot(gs[1])
pan1.imshow(odF2m[sortF2m][-1].transpose()[92:139,85:270],cmap='Blues',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
pan1.yaxis.set_ticks([16])  
pan1.set_yticklabels([])
pan1.xaxis.set_ticks([])        
pan1.spines["top"].set_visible(False)
pan1.spines["bottom"].set_visible(False) 
pan2=figure.add_subplot(gs[2])
pan2.imshow(odF2m[sortF2m][-1].transpose()[139:186,85:270],cmap='Greens',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
pan2.yaxis.set_ticks([18])  
pan2.set_yticklabels([])
pan2.xaxis.set_ticks([])        
pan2.spines["top"].set_visible(False)
pan2.spines["bottom"].set_visible(False) 
pan3=figure.add_subplot(gs[3])
pan3.imshow(odF2m[sortF2m][-1].transpose()[186:233,85:270],cmap='Reds',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
pan3.yaxis.set_ticks([20])
pan3.set_yticklabels([])  
pan3.xaxis.set_ticks([])        
pan3.spines["top"].set_visible(False)
pan3.spines["bottom"].set_visible(False) 
pan4=figure.add_subplot(gs[4])
pan4.imshow(odF2m[sortF2m][-1].transpose()[233:280,85:270],cmap='Magentas',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
pan4.yaxis.set_ticks([22])  
pan4.set_yticklabels([])
pan4.xaxis.set_ticks([])        
pan4.spines["top"].set_visible(False)

pan4.set_xticks([33,94,157])
pan4.tick_params(bottom=True,top=False)
pan4.set_xticklabels([-2,0,2])
pan4.set_xlabel(r'Wavevector $k_x$ [$k_{\rm{L}}$]')


figure.show()
plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure1v14c.pdf', dpi=500)
#plt.savefig('Z:/My Documents/papers etc/talks and posters/sample1.jpg', dpi=500)

#datafileF2m=np.load('28Feb2017_files_394-423.npz')
#
#odF2m=datafileF2m['odFiltAll']
#xCentersF2m=datafileF2m['xCenters']
#tlistF2m=datafileF2m['tlist']
#qlistF2m=datafileF2m['qlist']
#for i,q in enumerate(qlistF2m):
#    if q>-0.8:
#        qlistF2m[i] = qlistF2m[i] - 2.0
#        
#A,B,dA,dB = lineFit2.lineFit(tlistF2m,qlistF2m,'','')
#
#qlistF2mFit = tlistF2m*A
#sortF2m=np.argsort(qlistF2mFit)
#figure = plt.figure()
#figure.clear()
#figure.set_size_inches(3.0,1.7)
#
#gs = gridspec.GridSpec(5,1)
#
#gs.update(left=0.05, right=0.95, top=1.0, bottom = 0.2,wspace=0,hspace=0)
#
#vmax=0.4
#vmin=0.0
#
#pan0=figure.add_subplot(gs[0])
#pan0.imshow(odF2m[sortF2m][12].transpose()[45:92,85:270],cmap='Cyans',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
#pan0.yaxis.set_ticks([14]) 
#pan0.set_yticklabels([]) 
#pan0.xaxis.set_ticks([])        
#pan0.spines["bottom"].set_visible(False) 
#pan1=figure.add_subplot(gs[1])
#pan1.imshow(odF2m[sortF2m][12].transpose()[92:139,85:270],cmap='Blues',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
#pan1.yaxis.set_ticks([16])  
#pan1.set_yticklabels([])
#pan1.xaxis.set_ticks([])        
#pan1.spines["top"].set_visible(False)
#pan1.spines["bottom"].set_visible(False) 
#pan2=figure.add_subplot(gs[2])
#pan2.imshow(odF2m[sortF2m][12].transpose()[139:186,85:270],cmap='Greens',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
#pan2.yaxis.set_ticks([18])  
#pan2.set_yticklabels([])
#pan2.xaxis.set_ticks([])        
#pan2.spines["top"].set_visible(False)
#pan2.spines["bottom"].set_visible(False) 
#pan3=figure.add_subplot(gs[3])
#pan3.imshow(odF2m[sortF2m][12].transpose()[186:233,85:270],cmap='Reds',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
#pan3.yaxis.set_ticks([20])
#pan3.set_yticklabels([])  
#pan3.xaxis.set_ticks([])        
#pan3.spines["top"].set_visible(False)
#pan3.spines["bottom"].set_visible(False) 
#pan4=figure.add_subplot(gs[4])
#pan4.imshow(odF2m[sortF2m][12].transpose()[233:280,85:270],cmap='Magentas',vmin=vmin, vmax=vmax,interpolation='none',aspect =0.5)    
#pan4.yaxis.set_ticks([22])  
#pan4.set_yticklabels([])
#pan4.xaxis.set_ticks([])        
#pan4.spines["top"].set_visible(False)
#
#pan4.set_xticks([33,94,157])
#pan4.tick_params(bottom=True,top=False)
#pan4.set_xticklabels([-2,0,2])
#pan4.set_xlabel(r'Wavevector $k_x$ [$k_{\rm{L}}$]')
#figure.show()
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure1v13d.pdf', dpi=500)
#plt.savefig('Z:/My Documents/papers etc/talks and posters/sample2.jpg', dpi=500)