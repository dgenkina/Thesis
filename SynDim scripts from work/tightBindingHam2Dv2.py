# -*- coding: utf-8 -*-
"""
Created on Mon Apr 02 15:10:05 2018

@author: dng5
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 17:21:35 2018

@author: dng5
"""

import numpy as np
from numpy import linalg as LA
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#import readIgor

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
q=3.0
c=4.0/q#0.0#1064.0/790.0#
periodic=False
Flat=True

def TBham2d(tx,tm,p,q,kj,km):

    H=np.zeros((q,q), dtype=complex)
    kList = np.zeros(q)
    cent=np.int(q/2)
    for i in np.arange(q):
        kList[i]=kj-2.0*np.pi*p*(i-cent)/q
        H[i,i]=-2.0*tx*np.cos(kList[i])
        if i==0:
            if q==1:
                H[i,i]+=-2.0*tm*np.cos(km*q)
            else:
                H[i,i+1]=-tm
                H[i,q-1]=-tm*np.exp(1.0j*km*q)
        elif i==q-1:
            H[i,0]=-tm*np.exp(-1.0j*km*q)
            H[i,i-1]=-tm
        else:
            H[i,i+1]=-tm
            H[i,i-1]=-tm

    return H,kList
    
def getEigenspectrum(tx,tm,p,q, plot=True):
    kjList=np.linspace(-np.pi,np.pi,300)
    kmList=np.linspace(-np.pi/q,np.pi/q,300)
    Egrid=np.zeros((kjList.size,kmList.size,q))
    Vmag=np.zeros((kjList.size,kmList.size,q))
    mgrid=np.zeros((kjList.size,kmList.size))
    cent=np.int(q/2)
    mList=np.arange(q)-cent
    for kjind,kj in enumerate(kjList):
        for kmind, km in enumerate(kmList):
            H,kList=TBham2d(tx,tm,p,q,kj,km)
            E,V = sLA.eigh(H)
            Egrid[kjind,kmind,:] = E
            Vmag[kjind,kmind,:] = V[:,0]*np.conjugate(V[:,0])
            mgrid[kjind,kmind] = np.dot(mList,Vmag[kjind,kmind])

    if plot:        
        figure=plt.figure()
        pan=figure.add_subplot(111)
        pan.imshow(Egrid[:,:,0], cmap='Greys', extent=(kmList[0],kmList[-1],kjList[0],kjList[-1]))
    return Egrid,kjList,kmList,mgrid, Vmag
    

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


Egrid0,kj0,km0,mgrid0,Vmag0 = getEigenspectrum(0.5,0.5,0,1, plot=False)

EgridQ3P1,kjQ3P1,kmQ3P1,mgridQ3P1,VmagQ3P1 = getEigenspectrum(0.5,0.5,1,3, plot=False)
q3=3
cent=np.int(q3/2)
mList3=np.arange(q3)-cent
EgridQ5P2,kjQ5P2,kmQ5P2,mgridQ5P2,VmagQ5P2 = getEigenspectrum(0.5,0.5,2,5, plot=False)
q5=5
cent=np.int(q5/2)
mList5=np.arange(q5)-cent

figure = plt.figure()
figure.clear()
figure.set_size_inches(1.75,2.5)

gs = gridspec.GridSpec(3,1, height_ratios=[15, 5,3])

gs.update(left=0.2, right=0.98, top=0.95, bottom = 0.15, hspace=0.15,wspace=0.05)

s=km0.size
pan0=figure.add_subplot(gs[0])
pan0.imshow(Egrid0[:,:,0].transpose(), cmap='Greys',extent=(kj0[0]/(2.0*np.pi),kj0[-1]/(2.0*np.pi),km0[0]/(2.0*np.pi),km0[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_xticklabels([])
pan0.set_yticks([-0.5,-0.25,0.0,0.25,0.5])
pan0.set_yticklabels(['-1/2','-1/4','0','1/4','1/2'])


pan1=figure.add_subplot(gs[1])
pan1.imshow(EgridQ3P1[:,:,0].transpose(), cmap='Greys',extent=(kjQ3P1[0]/(2.0*np.pi),kjQ3P1[-1]/(2.0*np.pi),kmQ3P1[0]/(2.0*np.pi),kmQ3P1[-1]/(2.0*np.pi)))
pan1.set_xticks([-0.5,0.0,0.5])
pan1.set_xticklabels([])
pan1.set_yticks([-1.0/6.0,1.0/6.0])
pan1.set_yticklabels(['-1/6','1/6'])



pan2=figure.add_subplot(gs[2])
pan2.imshow(EgridQ5P2[:,:,0].transpose(), cmap='Greys',extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi),kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
pan2.set_xticks([-0.5,0.0,0.5])
pan2.set_yticks([-0.1,0.1])
pan2.set_yticklabels(['-1/10','1/10'])
pan2.set_xticklabels(['-1/2','0','1/2'])
pan2.set_xlabel(r'$q_x$ [$2k_L$]')


figure.text(0.05,0.6,r'$q_s$ [$2k_L$]', rotation='vertical', size=8)
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure3v6a.pdf', dpi=500)
#np.savez('figure3a',Egrid0=Egrid0[:,:,0].transpose(),kj0=kj0/(2.0*np.pi),km0=km0/(2.0*np.pi),
#         EgridQ3P1=EgridQ3P1[:,:,0].transpose(),kjQ3P1=kjQ3P1/(2.0*np.pi),kmQ3P1=kmQ3P1/(2.0*np.pi),
#         EgridQ5P2=EgridQ5P2[:,:,0].transpose(),kjQ5P2=kjQ5P2/(2.0*np.pi),kmQ5P2=kmQ5P2/(2.0*np.pi))


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

cmapDict={}
cmapDict[-2]='Magentas'
cmapDict[-1]='Reds'
cmapDict[0]='Greens'
cmapDict[1]='Blues'
cmapDict[2]='Cyans'
width=15
lw=10

figure = plt.figure()
figure.clear()
figure.set_size_inches(1.75,2.5)
gs = gridspec.GridSpec(2,1)

gs.update(left=0.2, right=0.95, top=0.95, bottom = 0.15, hspace=0.15,wspace=0.05)



pan1=figure.add_subplot(gs[0])
for ind in np.arange(q3):

    pan1.scatter(kjQ3P1/(2.0*np.pi),[ind-np.int(q3/2) for j in range(kmQ3P1.size)],s=width,c=VmagQ3P1[:,kmQ3P1.size/2,ind],vmin=0.0,vmax=1.0,cmap=cmapDict[ind-np.int(q3/2)], marker='_',linewidths=lw)
pan1.set_xticks([-0.5,0.0,0.5])
pan1.set_xticklabels([])
pan1.set_yticks([-1,0,1])
pan1.set_ylim(-2.5,2.5)
pan1.set_xlim(-0.5,0.5)



pan2=figure.add_subplot(gs[1])
for ind in np.arange(q5):

    pan2.scatter(kjQ5P2/(2.0*np.pi),[ind-np.int(q5/2) for j in range(kmQ5P2.size)],s=width,c=VmagQ5P2[:,kmQ5P2.size/2,ind],vmin=0.0,vmax=1.0,cmap=cmapDict[ind-np.int(q5/2)], marker='_',linewidths=lw)
pan2.set_xticks([-0.5,0.0,0.5])
pan2.set_yticks([-2,-1,0,1,2])
pan2.set_xticklabels(['-1/2','0','1/2'])
pan2.set_xlabel(r'$q_x$ [$2k_L$]')
pan2.set_ylim(-2.5,2.5)
pan2.set_xlim(-0.5,0.5)

figure.text(0.05,0.6,r'site $m$', rotation='vertical', size=8)
#plt.savefig('Z:/My Documents/papers etc/Bloch Osc/Figures/figure3v6b.pdf', dpi=500)

#np.savez('figure3b',kjQ3P1=kjQ3P1/(2.0*np.pi),VmagQ3P1=VmagQ3P1[:,kmQ3P1.size/2,:],
#         kjQ5P2=kjQ5P2/(2.0*np.pi),VmagQ5P2=VmagQ5P2[:,kmQ5P2.size/2,:])