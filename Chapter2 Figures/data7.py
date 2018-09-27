# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:47:49 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np
import ImagingModel3Fig7

tfinal = 0.0001
steps = 300
trange = np.linspace(0,tfinal,steps)
I0range = np.linspace(0.01,5.0,50)
nuS = np.zeros((np.size(I0range),steps))
nuM = np.zeros((np.size(I0range),steps))
nu1 = np.zeros((np.size(I0range),steps))
nu2 = np.zeros((np.size(I0range),steps))


for i in range(len(I0range)):
    inot=I0range[i]
    outputTuple = ImagingModel3Fig7.ImageS(inot,tfinal,steps,1.6)
    nuS[i,:]=outputTuple[1]
    outputTuple2 = ImagingModel3Fig7.ImageM(inot,tfinal,steps,1.6)
    nuM[i,:]=outputTuple2[1]
    nu1[i,:]=outputTuple2[2]
    nu2[i,:]=outputTuple2[3]
    
np.savez('fig7data.npz', I0range=I0range, nuS=nuS, nuM=nuM, nu1=nu1, nu2=nu2)