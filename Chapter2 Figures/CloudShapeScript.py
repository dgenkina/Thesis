# -*- coding: utf-8 -*-
"""
Created on Tue Jul 01 11:54:31 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np
import ImagingModel3Fig6 

"""This is a script to run the imaging model for a range of probe intensities"""

time = 0.0002 #s, so up to 200us
steps = 1000
I0range = np.linspace(0.1,5,40)
od0All = np.empty([I0range.size,steps])
od1All = np.empty([I0range.size,steps])
IfinalAll = np.empty([I0range.size,steps])
VatomAll = np.empty([I0range.size,steps+1])
widths = np.exp(np.linspace(np.log(0.0001e-8),np.log(1.0e-8),10))
ODofInotWidth = np.empty([I0range.size,widths.size])
zfinal = 0.0001
zrange = np.linspace(0,zfinal,100000)
gaussians = np.empty([zrange.size,widths.size])


def Simulate(width):
    for i in range(I0range.size):
        inot = I0range[i]
        outputTuple = ImagingModel3Fig6.Image(inot,time,steps,width)
    #    od=outputTuple[0]
        od0All[i,:] = outputTuple[1]
    return (od0All, outputTuple[5])
    #    od1All[i,:] = outputTuple[2]
    #    IfinalAll[i,:] = outputTuple[3]
    #    VatomAll[i,:] = outputTuple[4]
    
    #outfile = open("InotsScriptOutputs",'r+')
  #  np.savez("Data/SimulatedOD"+str(int(ODinit*1000)), od=od, od0All=od0All, od1All=od1All, IfinalAll=IfinalAll, VatomAll=VatomAll, I0range = I0range)
    #outfile.close()
    
for i in range(widths.size):
    width = widths[i]
    OutputTuple2 = Simulate(width)
    ODofInotWidth[:,i]=OutputTuple2[0][:,steps-1]
    gaussians[:,i]=OutputTuple2[1]
    np.savez("CloudWidths/Simulated2Clouds_"+str(i), gaussian=OutputTuple2[1], ODofInot=OutputTuple2[0], width=width)