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
figures in this paper https://arxiv.org/abs/1804.06345'''

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
grp2.attrs['Note']=r'''All data used in Figure2'''
grp2a=grp2.create_group('Figure2a')
grp2aFile=np.load('BandStructureFigureData.npz')
grp2a.attrs['Note']=r'''Theoretically calculated data for band structure figure 2a'''
grp2a['Energies']=grp2aFile['E']
grp2a['Energies'].attrs['Note']=r'''2-D array of energies with rows indexing different values of crystal momentum
and columns indexing different bands of the lattice'''
grp2a['crystalMomentumList']=grp2aFile['kList']/2.0
grp2a['crystalMomentumList'].attrs['Note']=r'''List of crystal momentum values used in figure, corresponding
to the rows of the Energy matrix'''
grp2a['modalPositions']=grp2aFile['modeTheory']
grp2a['modalPositions'].attrs['Note']=r'''2-D array of modal positions with rows indexing different values
of crystal momentum and columns indexing different bands of the lattice (same shape as Energies).
Represented by color in Figure2a'''
grp2b=grp2.create_group('Figure2b')
grp2b.attrs['Note']=r'''All absorption images, as 2-D arrays, used in figure 2b, along with the extracted
crystal momenta for each image'''
grp2bFile=np.load('dataForImagesFigure.npz')
grp2b['opticalDepth']=grp2bFile['odF2']
grp2b['crystalMomentumList']=grp2bFile['q']/2.0
grp2b['opticalDepth'].attrs['Note']=r'''A 3-D array of optical depths obtained
from absorption images presented in figure 2b. The first index indexes the image, 
in order of crystal momentum. The second and third indexes are the x and y axes of 
the image. In the figure, the images were sliced along the x (second index) direction
and stacked vertically to separate different spin components'''
grp2b['crystalMomentumList'].attrs['Note']=r'''Crystal momenta of each of the
images presented in Figure2b, corresponding to the first index of the 3-D optical
depth array'''
grp2c=grp2.create_group('Figure2c')
grp2c.attrs['Note']=r'''Data presented in all three panels of figure 2c'''
grp2cFile=np.load('magnetizationPlotData.npz')
grp2c0=grp2c.create_group('Figure2cTop')
grp2c0.attrs['Note']=r'''Data for top panel of figure 2c, chern number 0'''
grp2c0['crystalMomentumData']=grp2cFile['qlistData0']
grp2c0['crystalMomentumData'].attrs['Note']=r'''Crystal momentum for the experimental data points.''' 
grp2c0['crystalMomentumTheory']=grp2cFile['qlistTheory0']
grp2c0['crystalMomentumTheory'].attrs['Note']=r'''Crystal momentum for the red theory curve.''' 
grp2c0['modalPositionFromData']=grp2cFile['dataMode0']
grp2c0['modalPositionFromData'].attrs['Note']=r'''Modal position measured experimentally.''' 
grp2c0['ErrorInModalPositionFromData']=grp2cFile['dataModeError0']
grp2c0['ErrorInModalPositionFromData'].attrs['Note']=r'''Uncertainty in the modal position measured experimentally.''' 
grp2c0['modalPositionFromTheory']=grp2cFile['theoryMode0']
grp2c0['modalPositionFromTheory'].attrs['Note']=r'''Modal position calculated theoretically to give the red curve.''' 
grp2cP=grp2c.create_group('Figure2cMiddle')
grp2cP.attrs['Note']=r'''Data for middle panel of figure 2c, positive chern number'''
grp2cP['crystalMomentumData']=grp2cFile['qlistDataP']
grp2cP['crystalMomentumData'].attrs['Note']=r'''Crystal momentum for the experimental data points.''' 
grp2cP['crystalMomentumTheory']=grp2cFile['qlistTheoryP']
grp2cP['crystalMomentumTheory'].attrs['Note']=r'''Crystal momentum for the red theory curve.''' 
grp2cP['modalPositionFromData']=grp2cFile['dataModeP']
grp2cP['modalPositionFromData'].attrs['Note']=r'''Modal position measured experimentally.''' 
grp2cP['ErrorInModalPositionFromData']=grp2cFile['dataModeErrorP']
grp2cP['ErrorInModalPositionFromData'].attrs['Note']=r'''Uncertainty in the modal position measured experimentally.''' 
grp2cP['modalPositionFromTheory']=grp2cFile['theoryModeP']
grp2cP['modalPositionFromTheory'].attrs['Note']=r'''Modal position calculated theoretically to give the red curve.''' 
grp2cM=grp2c.create_group('Figure2cBottom')
grp2cM.attrs['Note']=r'''Data for middle panel of figure 2c, negative chern number'''
grp2cM['crystalMomentumData']=grp2cFile['qlistDataM']
grp2cM['crystalMomentumData'].attrs['Note']=r'''Crystal momentum for the experimental data points.''' 
grp2cM['crystalMomentumTheory']=grp2cFile['qlistTheoryM']
grp2cM['crystalMomentumTheory'].attrs['Note']=r'''Crystal momentum for the red theory curve.''' 
grp2cM['modalPositionFromData']=grp2cFile['dataModeM']
grp2cM['modalPositionFromData'].attrs['Note']=r'''Modal position measured experimentally.''' 
grp2cM['ErrorInModalPositionFromData']=grp2cFile['dataModeErrorM']
grp2cM['ErrorInModalPositionFromData'].attrs['Note']=r'''Uncertainty in the modal position measured experimentally.''' 
grp2cM['modalPositionFromTheory']=grp2cFile['theoryModeM']
grp2cM['modalPositionFromTheory'].attrs['Note']=r'''Modal position calculated theoretically to give the red curve.''' 


grp3=f.create_group('Figure3')
grp3a=f.create_group('Figure3a')
grp3a.attrs['Note']=r'''Theoretical calculation for the three panels of figure 3a'''
grp3aFile=np.load('figure3a.npz')
grp3a['EnergyTop']=grp3aFile['Egrid0']
grp3a['EnergyTop'].attrs['Note']=r'''For top plot, Lowest band energy as a function of crystal 
momentum q_x and crystal momentum q_s as a 2-D array, with q_x labelling the columns and q_s labeling the rows'''
grp3a['CrystalMomentumXtop']=grp3aFile['kj0']
grp3a['CrystalMomentumXtop'].attrs['Note']=r'''For top plot, List of crystal momentum values q_x'''
grp3a['CrystalMomentumStop']=grp3aFile['km0']
grp3a['CrystalMomentumStop'].attrs['Note']=r'''For top plot, List of crystal momentum values q_s'''
grp3a['EnergyMiddle']=grp3aFile['EgridQ3P1']
grp3a['EnergyMiddle'].attrs['Note']=r'''For middle plot, Lowest band energy as a function of crystal 
momentum q_x and crystal momentum q_s as a 2-D array, with q_x labelling the columns and q_s labeling the rows'''
grp3a['CrystalMomentumXMiddle']=grp3aFile['kjQ3P1']
grp3a['CrystalMomentumXMiddle'].attrs['Note']=r'''For middle plot, List of crystal momentum values q_x'''
grp3a['CrystalMomentumSMiddle']=grp3aFile['kmQ3P1']
grp3a['CrystalMomentumSMiddle'].attrs['Note']=r'''For middle plot, List of crystal momentum values q_s'''
grp3a['EnergyBottom']=grp3aFile['EgridQ5P2']
grp3a['EnergyBottom'].attrs['Note']=r'''For bottom plot, Lowest band energy as a function of crystal 
momentum q_x and crystal momentum q_s as a 2-D array, with q_x labelling the columns and q_s labeling the rows'''
grp3a['CrystalMomentumXBottom']=grp3aFile['kjQ5P2']
grp3a['CrystalMomentumXBottom'].attrs['Note']=r'''For bottom plot, List of crystal momentum values q_x'''
grp3a['CrystalMomentumSBottom']=grp3aFile['kmQ5P2']
grp3a['CrystalMomentumSBottom'].attrs['Note']=r'''For bottom plot, List of crystal momentum values q_s'''
grp3b=f.create_group('Figure3b')
grp3bFile=np.load('figure3b.npz')
grp3b.attrs['Note']=r'''Theoretical calculation for the two panels of figure 3b'''
grp3b['CrystalMomentumTop']=grp3bFile['kjQ3P1']
grp3b['CrystalMomentumTop'].attrs['Note']=r'''List of crystal momenta for top plot'''
grp3b['CrystalMomentumBottom']=grp3bFile['kjQ5P2']
grp3b['CrystalMomentumBottom'].attrs['Note']=r'''List of crystal momenta for bottom plot'''
grp3b['FractionalPopulationsTop']=grp3bFile['VmagQ3P1']
grp3b['FractionalPopulationsTop'].attrs['Note']=r'''2-D array containing fractional populations
in each m state for the top plot. The rows correspond to different crystal momenta and the columns
correspond to differen m states, from -1 to 1'''
grp3b['FractionalPopulationsBottom']=grp3bFile['VmagQ5P2']
grp3b['FractionalPopulationsBottom'].attrs['Note']=r'''2-D array containing fractional populations
in each m state for the Bottom plot. The rows correspond to different crystal momenta and the columns
correspond to differen m states, from -2 to 2'''

grp4=f.create_group('Figure4')
filenamelist=['07Mar2017_F1_chern_-1.npz','27Feb2017_F2_chern_-1.npz',
              '09Mar2017_Rf_Corrected.npz','22Mar2017_Rf_Corrected.npz',
              '08Mar2017_F1_chern_1.npz','28Feb2017_F2_chern_1.npz']
Slist=[1,2,1,2,1,2]
clist=[-1,-1,0,0,1,1]
for i,filename in enumerate(filenamelist):
    dataFile=np.load(filename)
    S=Slist[i]
    c=clist[i]
    stringShort=str(2*S+1)+'sitesChern'+str(c)
    string=str(2*S+1)+'-site ribbon, expected Chern number '+str(c)
    qlist=dataFile['qlist']
    fractionP=dataFile['fractionP']
    fraction0=dataFile['fraction0']
    fractionM=dataFile['fractionM']   
    fracs=np.array([fractionM,fraction0,fractionP])
    if S==2:
        fractionP2=dataFile['fractionP2']
        fractionM2=dataFile['fractionM2']
        fracs=np.array([fractionM2,fractionM,fraction0,fractionP,fractionP2])
    grp4['FractionalPopulations'+stringShort]=fracs
    grp4['FractionalPopulations'+stringShort].attrs['Note']=r'''2-D array of fractional populations,
    with rows indexing the m site (from lowest to highest) and columns indexing crystal momentum along x for the '''+ string
    grp4['CrystalMomentum'+stringShort]=fracs
    grp4['CrystalMomentum'+stringShort].attrs['Note']=r'''1-D array of crystal momenta along x,
    for the '''+ string
    

print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2)) 

f.close()