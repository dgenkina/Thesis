# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 15:38:31 2018

@author: swooty
"""

import numpy as np
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


#make the figure
plt.close(1)
fig = plt.figure(1,figsize=(6.2,6.0))
gs = gridspec.GridSpec(5,5)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.5,hspace=1.0)

#make ramp  on graph lines
tlist = np.linspace(0.0,25.0,1200)
latlist = np.zeros(tlist.size)
latlist[tlist<20.0] = tlist[tlist<20.0]*5.0/20.0
latlist[tlist>=20.0] = 5.0
latlist[tlist>23.0] = 0
OmegaList = np.zeros(tlist.size)
OmegaList[tlist>20.0] = (tlist[tlist>20.0]-20.0)*2.0/0.3
OmegaList[tlist>20.3] = 2
OmegaList[tlist>23] = 0
#make ramp on panel
panel=fig.add_subplot(gs[0,:])
panel.plot(tlist,latlist,'b-', label=r'$V_0$')
panel.plot(tlist,OmegaList,'r-',label = r'$\hbar\Omega$')
panel.text(21.5,2.8, r'$t$')
panel.set_xlim(0,25)
panel.arrow(x=20.3, y=2.5, dx=2.7,dy=0,length_includes_head=True,
            width=0.01,head_width=0.4,color='k')
panel.arrow(x=23, y=2.5, dx=-2.7,dy=0,length_includes_head=True,
            width=0.01,head_width=0.4,color='k')
panel.set_xlabel('Time [ms]')
panel.set_ylabel(r'Energy [$E_L$]')
panel.set_title('(a)')
plt.legend()


#make F=1 pulsing panel
panel = fig.add_subplot(gs[1:3,0:3])
fileF1 = np.load('SynDimPulsingF1.npz')
S=1
#plot data
panel.plot(fileF1['tList']*1.0e3,fileF1['fractionP'],'bo', label=r'$m_F$=+1')
panel.plot(fileF1['tList']*1.0e3,fileF1['fraction0'],'go', label=r'$m_F$=0')
panel.plot(fileF1['tList']*1.0e3,fileF1['fractionM'],'ro', label=r'$m_F$=-1')
#plot fit
panel.plot(fileF1['tForFit']*1.0e3,fileF1['pops_fitted'][S],'g-')
panel.plot(fileF1['tForFit']*1.0e3,fileF1['pops_fitted'][S+1],'b-')
panel.plot(fileF1['tForFit']*1.0e3,fileF1['pops_fitted'][S-1],'r-')
#panel.set_xlabel('Raman pulse time [ms]')
panel.set_ylabel('Fractional population')
plt.legend()
panel.set_title('(b)')
#print fitted values for omega and delta, uncertainties for caption
print fileF1['popt'], fileF1['pcov']


#make F=1 TOF image panel
panel = fig.add_subplot(gs[1:3,3:5])
filename = 'C:/Users/swooty/Documents/Thesis Data/Syn. dim. pulsing F=1/PIXIS_07Mar2017_0698.ibw'
od = readIgor.processIBW(filename)['rotOD'][460:660,385:630]
panel.imshow(od,vmin=-0.1,vmax=0.35)
panel.set_xticks([77,124,171])
panel.set_xticklabels([r'+1',r'0',r'-1'])
panel.axhline(y=100,xmin=0, xmax=0.5, color='white',linestyle='--')
panel.text(5,99,r'$0$',color='white')
panel.axhline(y=38,xmin=0, xmax=0.48, color='white',linestyle='--')
panel.text(5,37,r'$2k_L$',color='white')
panel.axhline(y=162,xmin=0, xmax=0.48, color='white',linestyle='--')
panel.text(5,161,r'$-2k_L$',color='white')
panel.axhline(y=18,xmin=0, xmax=0.72, color='white',linestyle='--')
panel.text(5,17,r'$2k_R$',color='white')
panel.axhline(y=186,xmin=0, xmax=0.3, color='white',linestyle='--')
panel.text(5,185,r'$-2k_R$',color='white')
panel.axhline(y=76,xmin=0.72, xmax=1.0, color='white',linestyle='--')
panel.text(150,75,r'$2k_R-2k_L$',color='white')
panel.axhline(y=124,xmin=0.3, xmax=1.0, color='white',linestyle='--')
panel.text(150,123,r'$2k_L-2k_R$',color='white')
panel.set_yticks([])
panel.set_title('(c)')
#panel.set_xlabel(r'position along ${\bf e}_s$')
panel.set_ylabel(r'momentum along ${\bf e}_x$')


#make F=2 pulsing panel
panel = fig.add_subplot(gs[3:5,0:3])
fileF2 = np.load('SynDimPulsingF2.npz')
S=2
#plot data
panel.plot(fileF2['tList']*1.0e3,fileF2['fractionP2'],'co', label=r'$m_F$=+2')
panel.plot(fileF2['tList']*1.0e3,fileF2['fractionP'],'bo', label=r'$m_F$=+1')
panel.plot(fileF2['tList']*1.0e3,fileF2['fraction0'],'go', label=r'$m_F$=0')
panel.plot(fileF2['tList']*1.0e3,fileF2['fractionM'],'ro', label=r'$m_F$=-1')
panel.plot(fileF2['tList']*1.0e3,fileF2['fractionM2'],'mo', label=r'$m_F$=-2')
#plot fit
panel.plot(fileF2['tForFit']*1.0e3,fileF2['pops_fitted'][S],'g-')
panel.plot(fileF2['tForFit']*1.0e3,fileF2['pops_fitted'][S+1],'b-')
panel.plot(fileF2['tForFit']*1.0e3,fileF2['pops_fitted'][S-1],'r-')
panel.plot(fileF2['tForFit']*1.0e3,fileF2['pops_fitted'][S+2],'c-')
panel.plot(fileF2['tForFit']*1.0e3,fileF2['pops_fitted'][S-2],'m-')
#print fitted values for omega and delta, uncertainties for caption
print fileF2['popt'], fileF2['pcov']
panel.set_xlabel(r'Raman pulse time $t$ [ms]')
panel.set_ylabel('Fractional population')
plt.legend(loc='right')
panel.set_title('(d)')

#make F=2 TOF image panel
panel = fig.add_subplot(gs[3:5,3:5])
filename = 'C:/Users/swooty/Documents/Thesis Data/Syn. dim. pulsing F=2/PIXIS_12Jul2017_0004.ibw'
od = readIgor.processIBW(filename)['rotOD'][430:630,360:605]
panel.imshow(od,vmin=-0.05,vmax=0.25)
panel.set_xticks([12,68,124,180,236])
panel.set_xticklabels([r'-2',r'-1',r'0',r'+1',r'+2'])
panel.set_yticks([])
panel.set_xlabel(r'position along ${\bf e}_s$')
panel.set_ylabel(r'momentum along ${\bf e}_x$')
panel.set_title('(e)')

plt.savefig('SynDimPulsing.pdf',transparent=True)