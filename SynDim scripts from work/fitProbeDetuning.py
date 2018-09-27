# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 14:21:28 2016

@author: dng5
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:20:22 2016

@author: dng5
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def lor(x,x0,Gamma,A,offset):
    L=A*(Gamma/2.0)**2.0/((x-x0)**2.0+(Gamma/2)**2.0) +offset
    return L 
roiList=np.array([[525, 555, 455, 485],
                  [520, 550, 510, 540],
                  [515, 545, 565, 595]])

fileroot='Y:/Data/2016/August/20/PIXIS_20Aug2016'  
key='probefreq'

#counts1, fractions1, waveDict1, probeAvg1=batchCountMultipleROIs(fileroot,64,113,roiList,background=True,bgndLoc="top")
counts2, fractions2, waveDict2, probeAvg2=batchCountMultipleROIs(fileroot,174,223,roiList,background=True,bgndLoc="top")
freqs=waveDict2[key]#np.append(waveDict1[key],waveDict2[key])
numP=counts2[:,0]#np.append(counts1[:,0],counts2[:,0])
numM=counts2[:,2]#np.append(counts1[:,2],counts2[:,2])

popt1,pcov1 = optimize.curve_fit(lor,freqs,numP, p0=(147.0,9.0,650.0,5.0))
print popt1, np.sqrt(np.diag(pcov1))
freqForFit=np.linspace(np.min(freqs),np.max(freqs),500)
data_fitted1=lor(freqForFit,*popt1)
popt2,pcov2 = optimize.curve_fit(lor,freqs,numM, p0=(147.0,9.0,650.0,5.0))
print popt2, np.sqrt(np.diag(pcov2))
data_fitted2=lor(freqForFit,*popt2)

fig=plt.figure()
pan=fig.add_subplot(1,1,1)
pan.plot(freqs,numP, 'bo')
pan.plot(freqForFit,data_fitted1,'b-')
pan.plot(freqs,numM, 'go')
pan.plot(freqForFit,data_fitted2,'g-')
pan.set_xlabel('Probe frequency [MHz]')
pan.set_ylabel('counted number')
pan.set_title('mF=+1 atoms (blue): x0=%.2f,Gamma=%.2f,A=%.1f,offset=%.1f,\n mF=-1 atoms(green): x0=%.2f,Gamma=%.2f,A=%.1f,offset=%.1f' %(popt1[0],popt1[1],popt1[2],popt1[3],popt2[0],popt2[1],popt2[2],popt2[3]))
