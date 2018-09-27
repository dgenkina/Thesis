# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 18:25:32 2014

@author: dng5
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import ImagingModel3Fig6


rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 10
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

width=1e-6
outputTuple = ImagingModel3Fig6.Image(1.2,0.0002,1000,width)
indeces = np.linspace(0,1000,6,dtype=int)
number=outputTuple[7]

figure1=plt.figure()
figure1.clear()
figure1.set_size_inches(5.5,3.5)
gs = gridspec.GridSpec(2,len(indeces)/2)
gs.update(left=0.2, right=0.95, top=0.85, bottom = 0.15)
figure1.suptitle(r'Cloud width =' +str(width*1e6)+ r' $\mu \rm{m}$' )
panel = {}
cmposition1 = []
cmvelocity1 = []
for i in range(len(indeces)):
    index = indeces[i]
    panel[i] = figure1.add_subplot(gs[i])
    pos=np.array(outputTuple[6])[:,index]
    avgpos=np.average(pos)
    cmposition1.append(avgpos)
    relpos=pos-avgpos
    vel=np.array(outputTuple[4])[:,index]
    avgvel=np.average(vel)
    cmvelocity1.append(avgvel)
    relvel=vel-avgvel
    panel[i].scatter(relpos*1e6,relvel, c=np.arange(number), cmap='jet',alpha=0.5)
    panel[i].set_xlim(-2,2)
    panel[i].set_ylim(-0.06,0.06)
    panel[i].set_title("t = "+str(index*0.2)+" us")
    panel[i].yaxis.set_ticks([])
    if i<3:
        panel[i].xaxis.set_ticks([])
    else:
        panel[i].xaxis.set_ticks([-2,-1,0,1,2])
panel[0].yaxis.set_ticks([-0.06,-0.03,0,0.03,0.06])
panel[3].yaxis.set_ticks([-0.06,-0.03,0,0.03,0.06])
panel[0].set_ylabel(r'$\Delta \rm{v}$ [m/s]')
panel[3].set_xlabel(r'$\Delta \rm{z}$ [$\mu \rm{m}$]')
figure1.show()
plt.savefig('figure6b.pdf', transparent=True)

width = 1e-5
outputTuple = ImagingModel3Fig6.Image(1.2,0.0002,1000,width)
indeces = np.linspace(0,1000,6,dtype=int)
number=outputTuple[7]

figure2=plt.figure()
figure2.clear()
figure2.set_size_inches(5.5,3.5)
gs = gridspec.GridSpec(2,len(indeces)/2)
gs.update(left=0.2, right=0.95, top=0.85, bottom = 0.15)
figure2.suptitle(r'Cloud width =' +str(width*1e6)+ r' $\mu \rm{m}$' )
panel = {}

cmposition2 = []
cmvelocity2 = []
for i in range(len(indeces)):
    index = indeces[i]
    panel[i] = figure2.add_subplot(gs[i])
    pos=np.array(outputTuple[6])[:,index]
    avgpos=np.average(pos)
    cmposition2.append(avgpos)
    relpos=pos-avgpos
    vel=np.array(outputTuple[4])[:,index]
    avgvel=np.average(vel)
    cmvelocity2.append(avgvel)
    relvel=vel-avgvel
    panel[i].scatter(relpos*1e6,relvel, c=np.arange(number),cmap='jet', alpha=0.5)
    panel[i].set_xlim(-20,20)
    panel[i].set_ylim(-0.4,0.4)
    panel[i].set_title("t = "+str(index*0.2)+" us")
    panel[i].yaxis.set_ticks([])
    if i<3:
        panel[i].xaxis.set_ticks([])
    else:
        panel[i].xaxis.set_ticks([-20,-10,0,10,20])
panel[0].yaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
panel[3].yaxis.set_ticks([-0.4,-0.2,0,0.2,0.4])
panel[0].set_ylabel(r'$\Delta \rm{v}$ [m/s]')
panel[3].set_xlabel(r'$\Delta \rm{z}$ [$\mu \rm{m}$]')
figure2.show()
plt.savefig('figure6a.pdf', transparent=True)


