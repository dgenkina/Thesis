# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 16:56:01 2017

@author: dng5
"""

#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np
import matplotlib.pyplot as plt


rcParams['axes.labelsize'] = 10
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8

rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

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

cdictC = {'red':   [(0.0,  0.0, 1.0), 
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)],

         'blue':  [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)]}
                   
cyans = mpl.colors.LinearSegmentedColormap('Cyans', cdictC)
plt.register_cmap(cmap=cyans)        

cdictM = {'red':   [(0.0,  0.0, 1.0), 
                   (1.0,  1.0, 0.0)],

         'green': [(0.0,  0.0, 1.0),
                   (1.0,  0.0, 1.0)],

         'blue':  [(0.0,  0.0, 1.0),
                   (1.0,  1.0, 1.0)]}
                   
magentas = mpl.colors.LinearSegmentedColormap('Magentas', cdictM)
plt.register_cmap(cmap=magentas)            
datafile2=np.load('SynDimBandStructure_F2_n7_Chern1.npz')
kList=datafile2['kList']
posInd=kList>0
negInd=kList<0

lw=10
width=10
S=2
fig = plt.figure()
fig.clear()
fig.set_size_inches(3.0,1.0)
pan=fig.add_subplot(111)      
pan.scatter(kList[posInd],[0-S for j in range(kList[posInd].size)],c=datafile2['pops'][posInd,0,0],cmap='Magentas',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[posInd],[1-S for j in range(kList[posInd].size)],c=datafile2['pops'][posInd,0,1],cmap='Reds',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[posInd],[2-S for j in range(kList[posInd].size)],c=datafile2['pops'][posInd,0,2],cmap='Greens',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[posInd],[3-S for j in range(kList[posInd].size)],c=datafile2['pops'][posInd,0,3],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[posInd],[4-S for j in range(kList[posInd].size)],c=datafile2['pops'][posInd,0,4],cmap='Cyans',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(datafile2['kList']+2.0,[0-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,0],cmap='Magentas',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(datafile2['kList']+2.0,[1-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,1],cmap='Reds',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(datafile2['kList']+2.0,[2-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,2],cmap='Greens',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(datafile2['kList']+2.0,[3-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,3],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(datafile2['kList']+2.0,[4-S for j in range(datafile2['kList'].size)],c=datafile2['pops'][:,0,4],cmap='Cyans',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[negInd]+4.0,[0-S for j in range(kList[negInd].size)],c=datafile2['pops'][negInd,0,0],cmap='Magentas',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[negInd]+4.0,[1-S for j in range(kList[negInd].size)],c=datafile2['pops'][negInd,0,1],cmap='Reds',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[negInd]+4.0,[2-S for j in range(kList[negInd].size)],c=datafile2['pops'][negInd,0,2],cmap='Greens',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[negInd]+4.0,[3-S for j in range(kList[negInd].size)],c=datafile2['pops'][negInd,0,3],cmap='Blues',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)
pan.scatter(kList[negInd]+4.0,[4-S for j in range(kList[negInd].size)],c=datafile2['pops'][negInd,0,4],cmap='Cyans',s=width,vmin=0.0,vmax=1.0/S,marker='_',linewidths=lw)

pan.set_xlabel(r'Crystal momentum $q_x$ [$k_L$]')

pan.yaxis.set_ticks([])

plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure1v3b.pdf', transparent=True)