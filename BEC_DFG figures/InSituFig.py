# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 15:33:10 2018

@author: swooty
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import optimize

from matplotlib import rc, rcParams
import readIgor

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


#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(4.5,6.0))
gs = gridspec.GridSpec(2,1)
gs.update(left=0.05, right=0.95, top=0.95, bottom = 0.05,wspace = 0.1,hspace=0.1)

def ThomasFermi(x,x0,omega,A,offset):
    tf = offset - A*(omega*(x-x0))**2.0
    out = np.zeros(x.size)
    out[tf > 0.0] = tf[tf > 0.0]
    out[tf <= 0.0] = 0.0
    return out

def Gaussian(x,x0,sigma,A,offset):
    return A*np.exp(-((x-x0)/sigma)**2.0) + offset
    
    
filename = 'C:/Users/swooty/Documents/Thesis Data/2017Aug20 Lattice calibrations/Flea3_20Aug2017_0074.ibw' 

outDict = readIgor.processIBW(filename)
od1 = outDict['od1'][330:380,250:340]
od2 = outDict['od2'][330:380,250:340]


#pan=figure.add_subplot(212)
#pan.imshow(od2)

linOD = np.sum(od1,axis=0)
xlist = np.arange(linOD.size)
TFtry = ThomasFermi(xlist,47.0,0.5,1.0,70.0)
Gtry = Gaussian(xlist,47.0,12.0,70.0,0.0)

popt,pcov = optimize.curve_fit(ThomasFermi, xlist,linOD, p0 = (47.0,0.5,1.0,70.0))
print popt,pcov
TFopt = ThomasFermi(xlist,*popt)

popt2,pcov2 = optimize.curve_fit(Gaussian, xlist,linOD, p0 = (47.0,12.0,70.0,0.0))
print popt2,pcov2
Gopt = Gaussian(xlist,*popt2)

pan=fig.add_subplot(gs[0])
pan.imshow(od1)
pan.set_xticks([])
pan.set_yticks([])
pan.set_title('(a)')

pan=fig.add_subplot(gs[1])
pan.plot(linOD,'b-')
pan.plot(xlist,TFopt,'r--')
pan.plot(xlist,Gopt, 'k-.')
pan.set_ylim(-3.0,72.0)
pan.set_xticks([])
pan.set_yticks([])
pan.set_title('(b)')
#pan=fig.add_subplot(212)
#pan.plot(np.sum(od2,axis=0),'b-')

plt.savefig('InSitu.pdf', transparent=True)