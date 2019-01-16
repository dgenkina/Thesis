# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 16:10:18 2019

@author: swooty
"""
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import optimize
import matplotlib as mpl


rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['axes.titlesize'] = 8


rcParams['text.usetex'] = True
rcParams['pdf.fonttype'] = 42 # True type fonts
rcParams['font.family'] = 'sans-serif'
rcParams['font.serif'] = ['']
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


rcParams['axes.linewidth'] = 0.75
rcParams['lines.linewidth'] = 0.75
rcParams['lines.markersize'] = 2

rcParams['xtick.major.size'] = 3      # major tick size in points
rcParams['xtick.minor.size'] = 2      # minor tick size in points
rcParams['xtick.major.width'] = 0.75       # major tick width in points
rcParams['xtick.minor.width'] = 0.75      # minor tick width in points

rcParams['ytick.major.size'] = 3      # major tick size in points
rcParams['ytick.minor.size'] = 2      # minor tick size in points
rcParams['ytick.major.width'] = 0.75       # major tick width in points
rcParams['ytick.minor.width'] = 0.75      # minor tick width in points

def gaussian(xlist,sigma,A,x0,B):
    return A*np.exp(-(xlist-x0)**2.0/(2.0*sigma**2.0))+B

def tangent(xlist,sigma,A,x0,B,x):
    slope = -A*(x-x0)*np.exp(-(x-x0)**2.0/(2.0*sigma**2.0))/(sigma**2.0)
    return slope*(xlist-x) + gaussian(x,sigma,A,x0,B)

xlist = np.linspace(-3.0,3.0,num=100)

ylist = gaussian(xlist,1.0,-1.0,0.0,0.0)

plt.close(1)
fig = plt.figure(1,figsize=(6.2,3.0))
gs = gridspec.GridSpec(1,3)
gs.update(left=0.05, right=0.95, top=0.85, bottom = 0.15,wspace = 0.1,hspace=0.2)
pan=fig.add_subplot(gs[0,0])
pan.plot(xlist,ylist,'b-')
pan.plot(0.0,gaussian(0.0,1.0,-1.0,0.0,0.0), 'ko', markersize=8 )
pan.set_xticks([])
pan.set_yticks([])
pan.set_title(r'(a) ${\bf F}=0$')

pan=fig.add_subplot(gs[0,1])
pan.plot(xlist,ylist,'b-')
xTanList = np.linspace(-2.0,0.0,num=100)
yTanList = tangent(xTanList,1.0,-1.0,0.0,0.0,-0.8)
pan.plot(-0.8,gaussian(-0.8,1.0,-1.0,0.0,0.0), 'ko', markersize=8 )
pan.plot(xTanList,yTanList,'k--')
pan.set_xticks([])
pan.set_yticks([])
pan.set_title(r'(b) ${\bf F}=|F_x|$')

pan=fig.add_subplot(gs[0,2])
pan.plot(xlist,ylist,'b-')
xTanList = np.linspace(0.0,2.0,num=100)
yTanList = tangent(xTanList,1.0,-1.0,0.0,0.0,0.8)
pan.plot(0.8,gaussian(0.8,1.0,-1.0,0.0,0.0), 'ko', markersize=8 )
pan.plot(xTanList,yTanList,'k--')
pan.set_title(r'(b) ${\bf F}=-|F_x|$')
pan.set_xticks([])
pan.set_yticks([])

plt.savefig('kickFigure.pdf',transparent=True)