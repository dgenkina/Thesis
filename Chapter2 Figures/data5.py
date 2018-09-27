# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 17:21:21 2014

@author: dng5
"""

import matplotlib.pyplot as plt
import numpy as np

tfinal = 0.0001
steps = 300
trange = np.linspace(0,tfinal,steps)
I0range = np.logspace(-2,1)
v = np.zeros((np.size(I0range),steps+1))



for i in range(len(I0range)):
    inot=I0range[i]
    outputTuple = Image(inot,tfinal,steps,1.6)
    v[i,:]=outputTuple[5][0]
    
np.savez('fig5data.npz', v=v, I0range=I0range)

