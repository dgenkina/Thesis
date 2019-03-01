# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:44:48 2019

@author: swooty
"""
import numpy as np
from scipy import linalg as sLA
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
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


cdict5 = {'red':   [(0.0,  1.0, 1.0),
                    (1.0/4.0,  1.0, 1.0),
                    (2.0/4.0,  0.0, 0.0),
                    (3.0/4.0,  0.0, 0.0),
                    (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
                    (1.0/4.0,  0.0, 0.0),
                    (2.0/4.0,  1.0, 1.0),
                    (3.0/4.0,  0.0, 0.0),
                    (1.0,  1.0, 1.0)],
                   
         'blue':  [(0.0,  1.0, 1.0),
                    (1.0/4.0,  0.0, 0.0),
                    (2.0/4.0,  0.0, 0.0),
                    (3.0/4.0,  1.0, 1.0),
                    (1.0,  1.0, 1.0)]}
fiveSite = mpl.colors.LinearSegmentedColormap('fiveSite', cdict5)
plt.register_cmap(cmap=fiveSite)      

cdict3 = {'red':   [(0.0,  1.0, 1.0),
                    (1.0/2.0,  0.0, 0.0),
                    (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
                    (1.0/2.0,  1.0, 1.0),
                    (1.0,  0.0, 0.0)],
                   
         'blue':  [(0.0,  0.0, 0.0),
                    (1.0/2.0,  0.0, 0.0),
                    (1.0,  1.0, 1.0)]}
threeSite = mpl.colors.LinearSegmentedColormap('threeSite', cdict5)
plt.register_cmap(cmap=threeSite)    
         
t=0.1

Egrid0,kj0,km0,mgrid0,Vmag0 = getEigenspectrum(t,t,0,1, plot=False)

EgridQ3P1,kjQ3P1,kmQ3P1,mgridQ3P1,VmagQ3P1 = getEigenspectrum(t,t,1,3, plot=False)
q3=3
cent=np.int(q3/2)
mList3=np.arange(q3)-cent
EgridQ5P2,kjQ5P2,kmQ5P2,mgridQ5P2,VmagQ5P2 = getEigenspectrum(t,t,2,5, plot=False)
q5=5
cent=np.int(q5/2)
mList5=np.arange(q5)-cent

figure = plt.figure(1)
figure.clear()
figure.set_size_inches(3.5,2.75)

gs = gridspec.GridSpec(1,2)
gs.update(left=0.4, right=0.95, top=0.95, bottom = 0.15, hspace=0.2,wspace=0.15)
 
pan0=figure.add_subplot(gs[0,0])
pan0.imshow(Egrid0[:,:,0].transpose(), cmap='Greys',extent=(kj0[0]/(2.0*np.pi),kj0[-1]/(2.0*np.pi),km0[0]/(2.0*np.pi),km0[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_yticks([-0.5,-0.25,0.0,0.25,0.5])
#pan0.set_xlabel(r'$q_x$ [$2k_L$]')
pan0.set_ylabel(r'$q_s$ [$2k_L$]')
pan0.set_yticklabels(['-1/2','-1/4','0','1/4','1/2'])
pan0.set_xticklabels([])

pan0=figure.add_subplot(gs[0,1])
pan0.imshow(mgrid0.transpose(), cmap='Greys',extent=(kj0[0]/(2.0*np.pi),kj0[-1]/(2.0*np.pi),km0[0]/(2.0*np.pi),km0[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_yticks([-0.5,-0.25,0.0,0.25,0.5])
pan0.set_yticklabels([])
pan0.set_xticklabels([])
#pan0.set_xlabel(r'$q_x$ [$2k_L$]')
#pan0.set_ylabel(r'$q_s$ [$2k_L$]')
fileroot = 'C:/Users/swooty/Documents/Work/Bloch Osc Paper/Figures/'
plt.savefig(fileroot+'TbBsMagP0Q1.pdf',transparent=True)

#figure = plt.figure(2)
#figure.clear()
#gs = gridspec.GridSpec(1,1)
#pan = figure.add_subplot(gs[0,0])
#pan.scatter(kj0,Egrid0[:,kj0.size/2,0],c=mgrid0[:,kj0.size/2],vmin=-1.0,vmax=1.0,marker='o')


figure = plt.figure(3)
figure.clear()
figure.set_size_inches(3.5,1.5)

gs = gridspec.GridSpec(1,2)
gs.update(left=0.4, right=0.95, top=0.95, bottom = 0.15, hspace=0.1,wspace=0.15)
 
pan0=figure.add_subplot(gs[0,0])
pan0.imshow(EgridQ3P1[:,:,0].transpose(), cmap='Greys',extent=(kjQ3P1[0]/(2.0*np.pi),kjQ3P1[-1]/(2.0*np.pi),kmQ3P1[0]/(2.0*np.pi),kmQ3P1[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_yticks([-1.0/6.0,1.0/6.0])
pan0.set_yticklabels(['-1/6','1/6'])
#pan0.set_xlabel(r'$q_x$ [$2k_L$]')
pan0.set_ylabel(r'$q_s$ [$2k_L$]')
pan0.set_xticklabels([])

pan0=figure.add_subplot(gs[0,1])
pan0.imshow(mgridQ3P1.transpose(), cmap='threeSite',extent=(kjQ3P1[0]/(2.0*np.pi),kjQ3P1[-1]/(2.0*np.pi),kmQ3P1[0]/(2.0*np.pi),kmQ3P1[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_yticks([-1.0/6.0,1.0/6.0])
pan0.set_yticklabels([])
pan0.set_xticklabels([])
#pan0.set_xlabel(r'$q_x$ [$2k_L$]')
#pan0.set_ylabel(r'$q_s$ [$2k_L$]')
plt.savefig(fileroot+'TbBsMagP1Q3.pdf',transparent=True)

#figure = plt.figure(4)
#figure.clear()
#gs = gridspec.GridSpec(1,1)
#pan = figure.add_subplot(gs[0,0])
#pan.scatter(kjQ3P1,EgridQ3P1[:,kjQ3P1.size/2,0],c=mgridQ3P1[:,kjQ3P1.size/2],vmin=-1.0,vmax=1.0,marker='o')

figure = plt.figure(5)
figure.clear()
figure.set_size_inches(3.5,1.0)

gs = gridspec.GridSpec(1,2)
gs.update(left=0.4, right=0.95, top=0.95, bottom = 0.15, hspace=0.3,wspace=0.15)
 
pan0=figure.add_subplot(gs[0,0])
pan0.imshow(EgridQ5P2[:,:,0].transpose(), cmap='Greys',extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi),kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_xticklabels(['-1/2','0','1/2'])
pan0.set_xlabel(r'$q_x$ [$2k_L$]')
pan0.set_ylabel(r'$q_s$ [$2k_L$]')
pan0.set_yticks([-0.1,0.1])
pan0.set_yticklabels(['-1/10','1/10'])

pan0=figure.add_subplot(gs[0,1])
pan0.imshow(mgridQ5P2.transpose(), cmap='fiveSite',extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi),kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
pan0.set_xticks([-0.5,0.0,0.5])
pan0.set_yticks([-0.1,0.1])
pan0.set_yticklabels([])
pan0.set_xticklabels(['-1/2','0','1/2'])
pan0.set_xlabel(r'$q_x$ [$2k_L$]')
#pan0.set_ylabel(r'$q_s$ [$2k_L$]')
plt.savefig(fileroot+'TbBsMagP2Q5.pdf',transparent=True)

#figure = plt.figure(6)
#figure.clear()
#gs = gridspec.GridSpec(1,1)
#pan = figure.add_subplot(gs[0,0])
#pan.scatter(kjQ5P2,EgridQ5P2[:,kjQ5P2.size/2,0],c=mgridQ5P2[:,kjQ5P2.size/2],vmin=-1.0,vmax=1.0,marker='o')


#s=5
#figure = plt.figure(3)
#figure.clear()
#figure.set_size_inches(2.5,2.5/np.float(s))
#gs = gridspec.GridSpec(2,1)
#gs.update(left=0.2, right=0.98, top=0.95, bottom = 0.15, hspace=0.0,wspace=0.0)
#
#kjmin = kjQ5P2[0]/(2.0*np.pi)
#kjmax = kjQ5P2[-1]/(2.0*np.pi)
#kjdel = (kjmax-kjmin)/np.float(s)
#N = kjQ5P2.size
#for i in range(s):
#pan1=figure.add_subplot(gs[0,:])
#pan1.imshow(EgridQ5P2[:,:,0].transpose(),cmap='Greys',extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi)/s,kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
#    print 'average:'
#    print str(np.average(mgridQ5P2.transpose()[:,i*N/s:(i+1)*N/s]))+' +/- '+str(np.std(mgridQ5P2.transpose()[:,i*N/s:(i+1)*N/s]))
#    print 'center:'
#    print mgridQ5P2.transpose()[N/2,i*N/s+N/s/2]
#    print 'm:'
#    m = np.int(np.round(np.average(mgridQ5P2.transpose()[:,i*N/s:(i+1)*N/s]),0))
#    print m
#rgba = fiveSite((mgridQ5P2.transpose()-np.min(mgridQ5P2))/(np.max(mgridQ5P2)-np.min(mgridQ5P2)))
#rgb = rgb=rgba[:,:,:3]
#hsv = mpl.colors.rgb_to_hsv(rgb)
#hsv[:,:,1] = (EgridQ5P2[:,:,0].transpose()-EgridQ5P2[:,:,0].min())/(EgridQ5P2[:,:,0].transpose().max()-EgridQ5P2[:,:,0].min())
#rgbNew = mpl.colors.hsv_to_rgb(hsv)    
#hsl = np.zeros(mgridQ5P2.transpose().shape + (3,))
#hsl[:, :, 1] = (EgridQ5P2[:,:,0].transpose()-EgridQ5P2[:,:,0].min())/(EgridQ5P2[:,:,0].transpose().max()-EgridQ5P2[:,:,0].min())
#hsl[:, :, 0] = (mgridQ5P2.transpose()+2.0)/4.0
#hsl[:, :, 2] = 0.5
#rgb = mpl.colors.hsv_to_rgb(hsl)
#d=pan1.imshow(rgbNew,extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi),kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))

#    pan1.imshow(EgridQ5P2[:,:,0].transpose()[:,i*N/s:(i+1)*N/s], cmap=cmapDict[m],extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi)/s,kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
#    pan1.set_xticks([])
#    if i==0:
#        pan1.set_yticks([-0.1,0.1])
#        pan1.set_yticklabels(['-1/10','1/10'])
#    else:    
#        pan1.set_yticks([])
##s=5
#kjmin = kjQ5P2[0]/(2.0*np.pi)
#kjmax = kjQ5P2[-1]/(2.0*np.pi)
#kjdel = (kjmax-kjmin)/np.float(s)
#N = kjQ5P2.size
#pan2=figure.add_subplot(gs[2])
#for i in range(s):
#    m = np.int(np.round(np.average(mgridQ3P1.transpose()[:,i*N/s:(i+1)*N/s]),0))
#    pan1.imshow(mgridQ5P2.transpose()[:,i*N/s:(i+1)*N/s], cmap=cmapDict[m],extent=(kjQ5P2[0]/(2.0*np.pi),kjQ5P2[-1]/(2.0*np.pi),kmQ5P2[0]/(2.0*np.pi),kmQ5P2[-1]/(2.0*np.pi)))
#pan2.set_xticks([-0.5,0.0,0.5])
#pan2.set_yticks([-0.1,0.1])
#pan2.set_yticklabels(['-1/10','1/10'])
#pan2.set_xticklabels(['-1/2','0','1/2'])
#pan2.set_xlabel(r'$q_x$ [$2k_L$]')
#
#
#figure.text(0.05,0.6,r'$q_s$ [$2k_L$]', rotation='vertical', size=8)


#
#s=3
#figure = plt.figure(2)
#figure.clear()
#figure.set_size_inches(2.5/np.float(s),2.5)
#gs = gridspec.GridSpec(1,s)
#gs.update(left=0.2, right=0.98, top=0.95, bottom = 0.15, hspace=0.0,wspace=0.0)
#
#kjmin = kjQ3P1[0]/(2.0*np.pi)
#kjmax = kjQ3P1[-1]/(2.0*np.pi)
#kjdel = (kjmax-kjmin)/np.float(s)
#N = kjQ3P1.size
#for i in range(s):
#    pan1=figure.add_subplot(gs[0,i])
#    print 'average:'
#    print str(np.average(mgridQ3P1.transpose()[:,i*N/s:(i+1)*N/s]))+' +/- '+str(np.std(mgridQ3P1.transpose()[:,i*N/s:(i+1)*N/s]))
#    print 'center:'
#    print mgridQ3P1.transpose()[N/2,i*N/s+N/s/2]
#    print 'm:'
#    m = np.int(np.round(np.average(mgridQ3P1.transpose()[:,i*N/s:(i+1)*N/s]),0))
#    print m
#    pan1.imshow(EgridQ3P1[:,:,0].transpose()[:,i*N/s:(i+1)*N/s], cmap=cmapDict[m],extent=(kjQ3P1[0]/(2.0*np.pi),kjQ3P1[-1]/(2.0*np.pi)/s,kmQ3P1[0]/(2.0*np.pi),kmQ3P1[-1]/(2.0*np.pi)))
#    pan1.set_xticks([])
#    if i==0:
#        pan1.set_yticks([-1.0/6.0,1.0/6.0])
#        pan1.set_yticklabels(['-1/6','1/6'])
#    else:    
#        pan1.set_yticks([])
#pan1.set_xticks([-0.5,0.0,0.5])
#pan1.set_xticklabels([])
#pan1.set_yticks([-1.0/6.0,1.0/6.0])
#pan1.set_yticklabels(['-1/6','1/6'])