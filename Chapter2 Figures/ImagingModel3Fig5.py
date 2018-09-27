# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:40:39 2014

@author: dng5
"""

#from pylab import *
import numpy as np
from numpy import pi

"""first define constants"""
hbar = 1.05457173*1e-34 #m^2kg/s
c = 299792485 #m/s

"""for K atoms, on resonant light, D2 line:"""
Gamma = 2*pi*6035000.0 #Hz
k = 2*pi*1304289.97000#1298518.93857 #m^-1 
omega = k*c
vrecoil = 0.01296541083 #m/s
isat = 17.5 #W/m^2
A = hbar*omega*Gamma/2/isat
B = k*vrecoil

"""Define superatom class"""
class SuperAtom:
    def __init__(self,position,velocity):
        self.z = position
        self.v = velocity
        self.zlist = [position]
        self.vlist = [velocity]
    def updateVelocity(self, Vadded):
        self.v += Vadded
        self.vlist.append(self.v)
    def updatePosition(self,time):
        self.z += self.v*time
        self.zlist.append(self.z)
        
        

def Image(Inot, tfinal,steps,ODinit):         
    """define initial atom distribution, and convert to superatoms.
    Assign initial position and velocity to each superatom"""
    zfinal = 0.0001 #m, so ~100um
    def n(x):
        return ODinit/zfinal/A
    
    zrange = np.linspace(0,zfinal,100000)
    dz = zrange[1]-zrange[0]#zfinal/100000.00
    
    superSize = 1e10 #number of atoms per unit area in one superatom
    integral = np.cumsum([n(z)*dz for z in zrange])/superSize
    od = ODinit #integral[zrange.size-1]*superSize*A
    superAtomNumber = int(integral[zrange.size-1]) 
    atomIndexes = np.arange(superAtomNumber)
    positions = np.interp(atomIndexes+1,integral,zrange)
    atoms = [SuperAtom(zed/2,0.0) for zed in positions]
    
    "define time grid and initialize intensity array"
   
    trange = np.linspace(0,tfinal,steps)
    dt = trange[1]-trange[0]#tfinal/float(steps)
 #   Inot = 1.3
    I = np.zeros([trange.size, superAtomNumber+1])
    I[:,0] = Inot
    
    for t in range(trange.size):
        sortedAtoms = atoms #sorted(atoms, key=lambda SuperAtom: SuperAtom.z)
        for atomIndex in atomIndexes:
            Delta = 2*k*sortedAtoms[atomIndex].v/Gamma #B*trange[t]*(I[t,atomIndex])/(1+I[t,atomIndex])
            rate = (I[t,atomIndex])/(1+(Delta**2)+I[t,atomIndex])
            I[t,atomIndex+1]=I[t,atomIndex]-A*superSize*rate
            sortedAtoms[atomIndex].updateVelocity(B*rate*dt*Gamma/2/k)
            sortedAtoms[atomIndex].updatePosition(dt)
        
    Ifinaltot = np.cumsum(I[:,superAtomNumber]*dt)/(trange+dt)
    nu = -np.log(Ifinaltot/Inot)     
    nu1 = od - (Inot-Ifinaltot)
    od2 = nu1 - ((B*(trange+dt))**2.0)*(1.0/(Ifinaltot+1.0)-1.0/(Inot+1.0)+np.log((Ifinaltot+1.0)/(Inot+1.0)))/3.0
    
    return (od, nu, nu1, od2, Ifinaltot, [atoms[i].vlist for i in range(superAtomNumber)], dz)



