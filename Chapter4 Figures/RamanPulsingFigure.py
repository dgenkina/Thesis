# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:06:25 2018

@author: swooty
"""

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
gs = gridspec.GridSpec(2,1)
gs.update(left=0.15, right=0.85, top=0.85, bottom = 0.15,wspace = 0.1)


#
#TOF image of F=1 Raman pulsed
filename = 'C:/Users/swooty/Documents/Thesis Data/Raman pulsing F=1/PIXIS_17Feb2016_0040.ibw'
od = readIgor.processIBW(filename)['rotOD'][410:680,475:765]
pan = fig.add_subplot(gs[0,0])
pan.imshow(od,vmin=-0.1, vmax=1.0)
pan.set_xticks([])
pan.set_ylabel(r'$k_x$')
pan.set_yticks([30,140,235])
pan.set_yticklabels([r'$2k_{\rm R}$',r'$0$',r'$-2k_{\rm R} $'])
pan.set_title('(a)')

pan.text(220,16,r'$-1$',color='white')
pan.text(143,16,r'$0$',color='white')
pan.text(20,19,r'$m_F=1$',color='white')

pan.text(10,250, 'SG gradient', color = 'white')
pan.arrow(10,260,260,0,color='white',width=0.1,head_width=10.0,length_includes_head = True,shape='full')

# Pulsing panel for F=1 Raman.
RamanPulsing = np.load('RamanPulsingF1.npz')

#Get fitted omega and delta in terms of raman recoils
print 'Omega:'
print RamanPulsing['popt'][0]
print '+/-'
print RamanPulsing['pcov'][0]

print 'delta:'
print RamanPulsing['popt'][1]
print '+/-'
print RamanPulsing['pcov'][1]

print 'epslion'
print 0.0207*ErecoilRaman/Erecoil

#make the F=1 pulsing panel
panel2 = fig.add_subplot(gs[1])

panel2.plot(RamanPulsing['time']*1.0e6,RamanPulsing['fractions'][:,0],'bo', label=r'$m_F$=+1')
panel2.plot(RamanPulsing['time']*1.0e6,RamanPulsing['fractions'][:,1],'go', label=r'$m_F$=0')
panel2.plot(RamanPulsing['time']*1.0e6,RamanPulsing['fractions'][:,2],'ro', label=r'$m_F$=-1')

panel2.plot(RamanPulsing['tForFit']*1.0e6,RamanPulsing['pops_fitted'][:,1],'g-')
panel2.plot(RamanPulsing['tForFit']*1.0e6,RamanPulsing['pops_fitted'][:,0],'b-')
panel2.plot(RamanPulsing['tForFit']*1.0e6,RamanPulsing['pops_fitted'][:,2],'r-')

panel2.set_xlabel(r'Raman pulse time [$\mu$s]')
panel2.set_ylabel('Fractional population')
panel2.set_title('(b)')
plt.legend()
fig.show()

plt.savefig('RamanPulsing.pdf', transparent=True)