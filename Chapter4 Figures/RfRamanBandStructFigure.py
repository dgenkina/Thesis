# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 19:22:47 2018

@author: swooty
"""

import matplotlib.pyplot as plt
#import numpy as np
import fitRaman
import fitRf



#import h5scripting
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np

# Use me to produce pdf files
# from matplotlib.backends.backend_pgf import FigureCanvasPgf
# matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

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
epsilon =0.04

#kList = np.linspace(-1.0,1.0,600)
#E00,m00 = fitRf.RfBandStructure(0.0, 0.0, epsilon, kList = kList, plot =False)
#E50,m50 = fitRf.RfBandStructure(5.0, 0.0, epsilon, kList = kList, plot =False)
#E01,m01 = fitRf.RfBandStructure(0.0, 1.0, epsilon, kList = kList, plot =False)
#E51,m51 = fitRf.RfBandStructure(5.0, 1.0, epsilon, kList = kList, plot =False)
#
#fig1 = plt.figure(1,figsize=(6.2,4.0))
#n=6
#gs = gridspec.GridSpec(2,2*n+1)
#gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace=2.0)
#
#
#pan1 = fig1.add_subplot(gs[0:n])
#for i in range(3):
#    d=pan1.scatter(kList,E00[:,i],c=m00[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
#pan1.set_xticks([])
#pan1.set_ylabel(r'Energy [$E_R$]')
#pan1.set_title('(a)')
#
#
#    
#pan2 = fig1.add_subplot(gs[n:2*n])
#for i in range(3):
#    pan2.scatter(kList,E50[:,i],c=m50[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
#pan2.set_xticks([])
#pan2.set_title('(b)')   
#    
#pan3 = fig1.add_subplot(gs[2*n+1:3*n+1])
#for i in range(3):
#    pan3.scatter(kList,E01[:,i],c=m01[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
#pan3.set_ylabel(r'Energy [$E_R$]')
#pan3.set_xlabel(r'$k_x$ [$k_R$]')
#pan3.set_title('(c)')  
#
#pan4 = fig1.add_subplot(gs[3*n+1:4*n+1])
#for i in range(3):
#    pan4.scatter(kList,E51[:,i],c=m51[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')    
#pan4.set_xlabel(r'$k_x$ [$k_R$]')
#pan4.set_title('(d)')  
#cbar_axes = fig1.add_subplot(gs[:,2*n])
#cbar = fig1.colorbar(d,cax=cbar_axes,ticks=np.array([-1,0,1]))
#
#plt.savefig('rfBandStructure.pdf', transparent=True)

kList = np.linspace(-3.0,3.0,1200)
E00,m00 = fitRaman.plotRamanBandStruct(0.0, 0.0, epsilon, kList = kList, plot =False)
E10,m10 = fitRaman.plotRamanBandStruct(1.0, 0.0, epsilon, kList = kList, plot =False)
E50,m50 = fitRaman.plotRamanBandStruct(5.0, 0.0, epsilon, kList = kList, plot =False)
E01,m01 = fitRaman.plotRamanBandStruct(0.0, 1.0, epsilon, kList = kList, plot =False)
E11,m11 = fitRaman.plotRamanBandStruct(1.0, 1.0, epsilon, kList = kList, plot =False)
E51,m51 = fitRaman.plotRamanBandStruct(5.0, 1.0, epsilon, kList = kList, plot =False)

fig1 = plt.figure(1,figsize=(6.2,4.0))
n=6
gs = gridspec.GridSpec(2,3*n+1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace=2.0)


pan1 = fig1.add_subplot(gs[0:n])
for i in range(3):
    d=pan1.scatter(kList,E00[:,i],c=m00[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
pan1.set_xticks([])
pan1.set_ylabel(r'Energy [$E_R$]')
pan1.set_title('(a)')
pan1.set_ylim([-3,12])

    
pan2 = fig1.add_subplot(gs[n:2*n])
for i in range(3):
    pan2.scatter(kList,E10[:,i],c=m10[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
pan2.set_xticks([])
pan2.set_title('(b)')
pan2.set_yticks([])
pan2.set_ylim([-3,12])   

pan3 = fig1.add_subplot(gs[2*n:3*n])
for i in range(3):
    pan3.scatter(kList,E50[:,i],c=m50[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
pan3.set_xticks([])
pan3.set_title('(c)')
pan3.set_yticks([])
pan3.set_ylim([-3,12])
    
pan4 = fig1.add_subplot(gs[3*n+1:4*n+1])
for i in range(3):
    pan4.scatter(kList,E01[:,i],c=m01[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')
pan4.set_ylabel(r'Energy [$E_R$]')
pan4.set_xlabel(r'$k_x$ [$k_R$]')
pan4.set_title('(d)') 
pan4.set_ylim([-3,12]) 

pan5 = fig1.add_subplot(gs[4*n+1:5*n+1])
for i in range(3):
    pan5.scatter(kList,E11[:,i],c=m11[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')    
pan5.set_xlabel(r'$k_x$ [$k_R$]')
pan5.set_title('(e)')
pan5.set_yticks([])  
pan5.set_ylim([-3,12])

pan6 = fig1.add_subplot(gs[5*n+1:6*n+1])
for i in range(3):
    pan6.scatter(kList,E51[:,i],c=m51[:,i],vmin=-1,vmax=1,cmap='jet_r',marker='.')    
pan6.set_xlabel(r'$k_x$ [$k_R$]')
pan6.set_title('(f)')  
pan6.set_yticks([])
pan6.set_ylim([-3,12])


cbar_axes = fig1.add_subplot(gs[:,3*n])
cbar = fig1.colorbar(d,cax=cbar_axes,ticks=np.array([-1,0,1]))

plt.savefig('RamanBandStructure.pdf', transparent=True)