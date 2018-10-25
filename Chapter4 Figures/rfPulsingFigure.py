# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 16:46:48 2018

@author: swooty
"""

import numpy as np
from scipy import optimize
from scipy import linalg as sLA
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib import rc, rcParams
import readIgor

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

hbar = 1.0545718e-34 # reduced Planck constant m^2 kg/s
mRb =1.44467e-25 #mass of rubidium in kg
lambdaR = 790e-9 # Raman wavelength in m
lambdaL = 1.064e-6 #lattice wavelength in m
Erecoil = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaL**2.0) #recoil energy
ErecoilRaman = (2.0*np.pi*hbar)**2.0/(2.0*mRb*lambdaR**2.0) 

#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,6.0))
gs = gridspec.GridSpec(3,2)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.1)


#
#TOF image of F=1 rf pulsed, 10Mar2017_0004-0033
filename = 'C:/Users/swooty/Documents/Thesis Data/10Mar2017 F=1 rf pulsing/PIXIS_10Mar2017_0004.ibw' 

od = readIgor.processIBW(filename,angle=-43)['rotOD'][515:595,370:630]
pan = fig.add_subplot(gs[0,0])
pan.imshow(od,vmin=-0.1, vmax=1.0)
pan.set_xticks([])
pan.set_yticks([])
pan.set_title('(a)')

pan.text(186,16,r'$-1$',color='white')
pan.text(143,16,r'$0$',color='white')
pan.text(48,19,r'$m_F=1$',color='white')

pan.text(10,60, 'SG gradient', color = 'white')
pan.arrow(10,70,240,0,color='white',width=0.1,head_width=10.0,length_includes_head = True,shape='full')

#Get a TOF image of F=2 rf pulse. 
filename = 'C:/Users/swooty/Documents/Thesis Data/2017Jan27 rf calibrations/PIXIS_27Jan2017_0035.ibw' 

od = readIgor.processIBW(filename,angle=-43)['rotOD'][515:595,370:630]
pan = fig.add_subplot(gs[0,1])
pan.imshow(od,vmin=-0.1, vmax=1.0)
pan.set_xticks([])
pan.set_yticks([])
pan.set_title('(b)')
pan.text(1,19,r'$m_F=-2$',color='white')
pan.text(82,16,r'$-1$',color='white')
pan.text(141,16,r'$0$',color='white')
pan.text(189,16,r'$1$',color='white')
pan.text(237,16,r'$2$',color='white')
pan.text(10,60, 'SG gradient', color = 'white')
pan.arrow(10,70,240,0,color='white',width=0.1,head_width=10.0,length_includes_head = True,shape='full')


# Pulsing panel for F=2 rf.
rfF2Pulsing = np.load('RfPulsingF2.npz')

#Get fitted omega and delta in terms of raman recoils
print 'Omega:'
print rfF2Pulsing['popt'][0]*ErecoilRaman/Erecoil 
print '+/-'
print  np.sqrt(np.diag(rfF2Pulsing['pcov']))[0]*ErecoilRaman/Erecoil

print 'delta:'
print rfF2Pulsing['popt'][1]*ErecoilRaman/Erecoil 
print '+/-'
print  np.sqrt(np.diag(rfF2Pulsing['pcov']))[1]*ErecoilRaman/Erecoil

print 'epslion'
print 0.0207*ErecoilRaman/Erecoil

# make the panel for F=2 pulsing data
S=2
panel = fig.add_subplot(gs[2,:])
#panel.set_title( r'$\Omega$ = ' + str(np.round(popt[0],2)) + r' $E_L/h$,$\delta$ = '+str(np.round(popt[1],3))+r' $E_L/h$')#, U= '+str(np.round(popt[2],3))+r'$E_L$')
panel.plot(rfF2Pulsing['tList']*1.0e3,rfF2Pulsing['fractionP2'],'co', label=r'$m_F$=+2')
panel.plot(rfF2Pulsing['tList']*1.0e3,rfF2Pulsing['fractionP'],'bo', label=r'$m_F$=+1')
panel.plot(rfF2Pulsing['tList']*1.0e3,rfF2Pulsing['fraction0'],'go', label=r'$m_F$=0')
panel.plot(rfF2Pulsing['tList']*1.0e3,rfF2Pulsing['fractionM'],'ro', label=r'$m_F$=-1')
panel.plot(rfF2Pulsing['tList']*1.0e3,rfF2Pulsing['fractionM2'],'mo', label=r'$m_F$=-2')


panel.plot(rfF2Pulsing['tForFit']*1.0e3,rfF2Pulsing['pops_fitted'][S],'g-')
panel.plot(rfF2Pulsing['tForFit']*1.0e3,rfF2Pulsing['pops_fitted'][S+1],'b-')
panel.plot(rfF2Pulsing['tForFit']*1.0e3,rfF2Pulsing['pops_fitted'][S-1],'r-')
panel.plot(rfF2Pulsing['tForFit']*1.0e3,rfF2Pulsing['pops_fitted'][S+2],'c-')
panel.plot(rfF2Pulsing['tForFit']*1.0e3,rfF2Pulsing['pops_fitted'][S-2],'m-')
panel.set_xlabel('Rf pulse time [ms]')
panel.set_ylabel('Fractional population')
panel.set_title('(d)')
plt.legend()

#Get F=1 pulsing data
rfF1Pulsing = np.load('RfPulsingF1.npz')
S=1

#Pring F=1 pulsing fitted omega and delta in terms of raman recoils
print 'Omega:'
print rfF1Pulsing['popt'][0]*ErecoilRaman/Erecoil 
print '+/-'
print  np.sqrt(np.diag(rfF1Pulsing['pcov']))[0]*ErecoilRaman/Erecoil

print 'delta:'
print rfF1Pulsing['popt'][1]*ErecoilRaman/Erecoil 
print '+/-'
print  np.sqrt(np.diag(rfF1Pulsing['pcov']))[1]*ErecoilRaman/Erecoil

#make the F=1 pulsing panel
panel2 = fig.add_subplot(gs[1,:])

panel2.plot(rfF1Pulsing['tList']*1.0e3,rfF1Pulsing['fractionP'],'bo', label=r'$m_F$=+1')
panel2.plot(rfF1Pulsing['tList']*1.0e3,rfF1Pulsing['fraction0'],'go', label=r'$m_F$=0')
panel2.plot(rfF1Pulsing['tList']*1.0e3,rfF1Pulsing['fractionM'],'ro', label=r'$m_F$=-1')

panel2.plot(rfF1Pulsing['tForFit']*1.0e3,rfF1Pulsing['pops_fitted'][S],'g-')
panel2.plot(rfF1Pulsing['tForFit']*1.0e3,rfF1Pulsing['pops_fitted'][S+1],'b-')
panel2.plot(rfF1Pulsing['tForFit']*1.0e3,rfF1Pulsing['pops_fitted'][S-1],'r-')

#panel2.set_xlabel('Rf pulse time [ms]')
panel2.set_ylabel('Fractional population')
panel2.set_title('(c)')
plt.legend()
fig.show()

plt.savefig('RfPulsing.pdf', transparent=True)