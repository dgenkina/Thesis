# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:25:01 2017

@author: dng5
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

h5name = "Measuring_Topology_in_Hofstadter_Ribbons_data.hdf5"

FileDocString='''This file contains all the data that was used to create 
figures in this paper '''

f=h5py.File(h5name,'a')
print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2))
#Clear file
for key in f.keys():
    f.__delitem__(key)
print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2))   


f.attrs['Note']=FileDocString


grp1=f.create_group('Figure1')
grp1.attrs['Note']=r'''Data displayed in Figure 1'''
grp1['Figure1c']=np.load('28Feb2017_files_454-483.npz')['odFiltAll'][20]
grp1['Figure1c'].attrs['Note']=r'''2D array containing uncropped time of flight image used in Figure 1c'''
grp1['Figure1d']=np.load('28Feb2017_files_394-423.npz')['odFiltAll'][13]
grp1['Figure1d'].attrs['Note']=r'''2D array containing uncropped time of flight image used in Figure 1d'''

grp2=f.create_group('Figure2')
grp2.attrs['Note']


datafilePulsing1=np.load('22Aug2017_files_4-33.npz')
GroupNamePulsing1='PulsingData1'
grpF2pulsing1=f.create_group(GroupNamePulsing1)
grpF2pulsing1.attrs['Note']=r'''Sample pulsing data from 22Aug2017, files 22Aug2017_4-34. 
The lattice depth here was 4.4 E_L and the 
fitted Raman coupling strenth was 0.62 E_L with a fitted detuning of -0.004 E_L. Quadratic zeeman shift was 0.02 E_L
The lattice was turned on adiabatically, and the RAman was turned on in 300 us'''
grpF2pulsing1['tlist']=datafilePulsing1['tlist']
grpF2pulsing1['tlist'].attrs['Note']=r'Raman pulse time in seconds'
grpF2pulsing1['fractionP2']=datafilePulsing1['fractionP2']
grpF2pulsing1['fractionP2'].attrs['Note']=r'fraction in mF=+2 state'
grpF2pulsing1['fractionP']=datafilePulsing1['fractionP']
grpF2pulsing1['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF2pulsing1['fraction0']=datafilePulsing1['fraction0']
grpF2pulsing1['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF2pulsing1['fractionM']=datafilePulsing1['fractionM']
grpF2pulsing1['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF2pulsing1['fractionM2']=datafilePulsing1['fractionM2']
grpF2pulsing1['fractionM2'].attrs['Note']=r'fraction in mF=-2 state'
grpF2pulsing1['kgrid']=np.fft.fftfreq(datafilePulsing1['tlist'].size,d=(datafilePulsing1['tlist'][1]-datafilePulsing1['tlist'][0]))
grpF2pulsing1['kgrid'].attrs['Note']=r'list of k states for fft'
grpF2pulsing1['fft0']=np.fft.fft(datafilePulsing1['fraction0'])
grpF2pulsing1['fft0'].attrs['Note']=r'fft of fractional population in mF=0 state'
grpF2pulsing1['fftP']=np.fft.fft(datafilePulsing1['fractionP'])
grpF2pulsing1['fftP'].attrs['Note']=r'fft of fractional population in mF=1 state'
grpF2pulsing1['fftM']=np.fft.fft(datafilePulsing1['fractionM'])
grpF2pulsing1['fftM'].attrs['Note']=r'fft of fractional population in mF=-1 state'
grpF2pulsing1['fftP2']=np.fft.fft(datafilePulsing1['fractionP2'])
grpF2pulsing1['fftP2'].attrs['Note']=r'fft of fractional population in mF=2 state'
grpF2pulsing1['fftM2']=np.fft.fft(datafilePulsing1['fractionM2'])
grpF2pulsing1['fftM2'].attrs['Note']=r'fft of fractional population in mF=-2 state'
grpF2pulsing1['Image']=np.load('22Aug2017_0006_processedImage.npz')['odFiltered']
grpF2pulsing1['Image'].attrs['Note']=r'''Od from file 23Aug2017_0006. Lattice orders are along the vertical axis,
and spin states are along horizontal axis, with mF=-2 being all the way on the left'''


print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2)) 

f.close()