# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 13:08:23 2014

@author: dng5
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import rc, rcParams
import numpy as np
import scipy
from scipy import special

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

ODbest = np.load('IsatCounts.npz')['ODbest']
IsatCounts = np.load('IsatCounts.npz')['IsatCounts']

simul = np.load("SimulatedOD"+str(int(ODbest*1000))+".npz")

times = (40.00, 75.00, 100.00)
probeCl={}
Rod0Clean={}
for time in times:
    probeCl[time]=np.load('CleanODProbe_'+str(time)+'.npz')['probeCl']
    Rod0Clean[time]=np.load('CleanODProbe_'+str(time)+'.npz')['Rod0Cl']
    
    
timeIndexes = (200, 375, 500)
odfigure = plt.figure()
odfigure.clear()
odfigure.set_size_inches(3.2,3.0)
gs = gridspec.GridSpec(1,1)
gs.update(left=0.15, right=0.95, top=0.95, bottom = 0.15, hspace=0.1, wspace=0.25)
odpanel = odfigure.add_subplot(gs[0])
color={}
color[0]='b'
color[1]='g'
color[2]='r'
for t in range(len(timeIndexes)):
    index = timeIndexes[t]
    odpanel.plot(probeCl[times[t]]/times[t]/IsatCounts,  Rod0Clean[times[t]] , color[t]+'o', ms=3)
    odpanel.plot(simul['I0range'], simul['od0All'][:,index], color[t]+'-',linewidth=2)
    odpanel.set_xlabel(r'$I_0/I_{\rm{sat}}$')
    odpanel.set_ylabel(r'$OD$')
#odpanel.text(2.6,0.38,'t=40us', rotation=-7, size=10 )
#odpanel.text(2.6,0.275,'t=75us', rotation=-5, size=10 )
#odpanel.text(2.6,0.13,'t=100us',rotation=-5, size=10 )
odpanel.set_xlim(0,3)
#odpanel.set_title('Isat = '+str(IsatCounts)+' counts/us')
odfigure.show()
plt.savefig('figure8.pdf', transparent=True)

"""ODbest = np.load('IsatCountsModel.npz')['ODbest']
IsatCounts = np.load('IsatCountsModel.npz')['IsatCounts']

odfigure = plt.figure()
odpanel = odfigure.add_subplot(1,1,1)
color={}
color[0]='b'
color[1]='g'
color[2]='r'
for t in range(len(timeIndexes)):
    index = timeIndexes[t]
    odpanel.plot(probeCl[times[t]]/times[t]/IsatCounts,  Rod0Clean[times[t]] , color[t]+'o')
    odpanel.plot(simul['I0range'], [ODbest +scipy.special.lambertw(inot*np.exp(inot-ODbest))-inot for inot in simul['I0range']], color[t]+'-',linewidth=2)
    odpanel.set_xlabel(r'$I_0/I_{sat}$', size=15)
    odpanel.set_ylabel(r'$\nu$', size=15)
odpanel.legend((odpanel.lines[0],odpanel.lines[2],odpanel.lines[4]),('t=40us','t=75us','t=100us'))
odpanel.set_xlim(0,3)
#odpanel.set_title('Isat = '+str(IsatCounts)+' counts/us')
odfigure.show()"""