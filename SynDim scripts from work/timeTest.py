# -*- coding: utf-8 -*-
"""
Created on Tue May 01 11:57:28 2018

@author: dng5
"""


import numpy as np
import time


d=100
A=np.random.rand(d,d)
u=np.random.rand(d)
i=0
print 'Initialized vectors'

t1=time.clock()
while i<1000:
    c=np.dot(A,u)
    i = (i
    +1)
t2=time.clock()
print 'For dimension d = '+str(d)+' 1000 matrix multiplications took ' + str(t2-t1) +' s'