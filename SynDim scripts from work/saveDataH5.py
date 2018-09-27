# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:08:51 2017

@author: dng5
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

h5name = "2017_08_25_BlochOscFiguresT2.hdf5"

FileDocString='Data for Bloch oscillations in a synthetic dimension project'

f=h5py.File(h5name,'a')
print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2))
#Clear file
for key in f.keys():
    f.__delitem__(key)
print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2))   


f.attrs['Note']=FileDocString

datafile=np.load('SynDimBandStructure_F2_n7_Chern1Flat.npz')
GroupName1='BandStructure_F2'
docstring1=r'''Data for band structure plot as well as populations in the lowest band.
Here, the Raman coupling is flat accross spin states, ie no CG coefficients,
F=2,number of lattice orders included is 15, Raman coupling strength 
$\Omega$=%.3f $E_L$, detuning $\delta$=%.3f, and quadratic shift $\epsilon$=%.3f'''%(datafile['omega'],datafile['delta'],datafile['epsilon'])
grp1 = f.create_group(GroupName1)
grp1.attrs['Note']=docstring1
grp1["kList"] = datafile['kList']
grp1["kList"].attrs['Note'] = r"List of quasimomenta used in units of $k_L$"
grp1["E"] = datafile['E'][:,:5]
grp1["E"].attrs['Note'] = r"Energies of the lowest 5 bands in units of $E_L$. Energy in ith band is E[:,i]"
grp1["m"] = datafile['m'][:,:5]
grp1["m"].attrs['Note'] = r"Magnetizations of the lowest 5 bands. Magnetization in ith band is m[:,i]"
grp1["pops"]= datafile['pops'][:,0,:]
grp1["pops"].attrs['Note']=r"Lowest band populations in each spin state. Population in each spin state mF is pops[:,mF-2]"
grp1.attrs['makeBandStructureFigure']= '''fig=plt.figure()
panel=fig.add_subplot(1,1,1)
for i in range(5):   
d=panel.scatter(grp1['kList'],grp1['E'][:,i],c=grp1['m'][:,i],vmin=-2,vmax=2, marker='_')
panel.set_xlabel(r'$q/k_L$')
panel.set_ylabel(r'$E/E_L$')
plt.colorbar(d)'''
grp1.attrs['makeLowestBandFigure']= '''fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)   
pan2.set_title('Lowest band')    
for i in range(5):    
    pan2.scatter(grp1['kList'],[i for j in range(grp1['kList'].size)],c=grp1['pops'][:,i],cmap='Blues',vmin=0.0,vmax=1.0,marker='_',linewidths=10)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'$q/k_L$')'''


datafileF2r=np.load('SynDimBandStructure_F2_n7_Chern1.npz')
GroupNameF2r='BandStructure_F2_real'
docstringF2r=r'''Data for band structure plot as well as populations in the lowest band.
Here, the Raman coupling is spin dependent, ie CG coefficients are accounted for,
F=2,number of lattice orders included is 15, Raman coupling strength 
$\Omega$=%.3f $E_L$, detuning $\delta$=%.3f, and quadratic shift $\epsilon$=%.3f'''%(datafileF2r['omega'],datafileF2r['delta'],datafileF2r['epsilon'])
grpF2r = f.create_group(GroupNameF2r)
grpF2r.attrs['Note']=docstringF2r
grpF2r["kList"] = datafileF2r['kList']
grpF2r["kList"].attrs['Note'] = r"List of quasimomenta used in units of $k_L$"
grpF2r["E"] = datafileF2r['E'][:,:5]
grpF2r["E"].attrs['Note'] = r"Energies of the lowest 5 bands in units of $E_L$. Energy in ith band is E[:,i]"
grpF2r["m"] = datafileF2r['m'][:,:5]
grpF2r["m"].attrs['Note'] = r"Magnetizations of the lowest 5 bands. Magnetization in ith band is m[:,i]"
grpF2r["pops"]= datafileF2r['pops'][:,0,:]
grpF2r["pops"].attrs['Note']=r"Lowest band populations in each spin state. Population in each spin state mF is pops[:,mF-2]"
grpF2r.attrs['makeBandStructureFigure']= '''fig=plt.figure()
panel=fig.add_subplot(1,1,1)
for i in range(5):   
d=panel.scatter(grpF2r['kList'],grpF2r['E'][:,i],c=grpF2r['m'][:,i],vmin=-2,vmax=2, marker='_')
panel.set_xlabel(r'$q/k_L$')
panel.set_ylabel(r'$E/E_L$')
plt.colorbar(d)'''
grpF2r.attrs['makeLowestBandFigure']= '''fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)   
pan2.set_title('Lowest band')    
for i in range(5):    
    pan2.scatter(grpF2r['kList'],[i for j in range(grpF2r['kList'].size)],c=grpF2r['pops'][:,i],cmap='Blues',vmin=0.0,vmax=1.0,marker='_',linewidths=10)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'$q/k_L$')'''

datafile2=np.load('SynDimBandStructure_F1_n7_Chern1Flat.npz')
GroupName2='BandStructure_F1'
docstring2=r'''Data for band structure plot as well as populations in the lowest band.
Here, the Raman coupling is flat accross spin states, ie no CG coefficients,
F=1,number of lattice orders included is 15, Raman coupling strength 
$\Omega$=%.3f $E_L$, detuning $\delta$=%.3f, and quadratic shift $\epsilon$=%.3f'''%(datafile2['omega'],datafile2['delta'],datafile2['epsilon'])
grp2 = f.create_group(GroupName2)
grp2.attrs['Note']=docstring2
grp2["kList"] = datafile2['kList']
grp2["kList"].attrs['Note'] = r"List of quasimomenta used in units of $k_L$"
grp2["E"] = datafile2['E'][:,:3]
grp2["E"].attrs['Note'] = r"Energies of the lowest 3 bands in units of $E_L$. Energy in ith band is E[:,i]"
grp2["m"] = datafile2['m'][:,:3]
grp2["m"].attrs['Note'] = r"Magnetizations of the lowest 3 bands. Magnetization in ith band is m[:,i]"
grp2["pops"]= datafile2['pops'][:,0,:]
grp2["pops"].attrs['Note']=r"Lowest band populations in each spin state. Population in each spin state mF is pops[:,mF-1]"
grp2.attrs['makeBandStructureFigure']= '''fig=plt.figure()
panel=fig.add_subplot(1,1,1)
for i in range(3):   
d=panel.scatter(grp2['kList'],grp2['E'][:,i],c=grp2['m'][:,i],vmin=-1,vmax=1, marker='_')
panel.set_xlabel(r'$q/k_L$')
panel.set_ylabel(r'$E/E_L$')
plt.colorbar(d)'''
grp2.attrs['makeLowestBandFigure']= '''fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)   
pan2.set_title('Lowest band')    
for i in range(3):    
    pan2.scatter(grp2['kList'],[i for j in range(grp2['kList'].size)],c=grp2['pops'][:,i],cmap='Blues',vmin=0.0,vmax=1.0,marker='_',linewidths=10)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'$q/k_L$')'''

datafileF1r=np.load('SynDimBandStructure_F1_n7_Chern1.npz')
GroupNameF1r='BandStructure_F1_real'
docstringF1r=r'''Data for band structure plot as well as populations in the lowest band.
Here, the Raman coupling is spin dependent, ie CG coefficients are accounted for,
F=1,number of lattice orders included is 15, Raman coupling strength 
$\Omega$=%.3f $E_L$, detuning $\delta$=%.3f, and quadratic shift $\epsilon$=%.3f'''%(datafileF1r['omega'],datafileF1r['delta'],datafileF1r['epsilon'])
grpF1r = f.create_group(GroupNameF1r)
grpF1r.attrs['Note']=docstringF1r
grpF1r["kList"] = datafileF1r['kList']
grpF1r["kList"].attrs['Note'] = r"List of quasimomenta used in units of $k_L$"
grpF1r["E"] = datafileF1r['E'][:,:3]
grpF1r["E"].attrs['Note'] = r"Energies of the lowest 3 bands in units of $E_L$. Energy in ith band is E[:,i]"
grpF1r["m"] = datafileF1r['m'][:,:3]
grpF1r["m"].attrs['Note'] = r"Magnetizations of the lowest 3 bands. Magnetization in ith band is m[:,i]"
grpF1r["pops"]= datafileF1r['pops'][:,0,:]
grpF1r["pops"].attrs['Note']=r"Lowest band populations in each spin state. Population in each spin state mF is pops[:,mF-1]"
grpF1r.attrs['makeBandStructureFigure']= '''fig=plt.figure()
panel=fig.add_subplot(1,1,1)
for i in range(3):   
d=panel.scatter(grpF1r['kList'],grpF1r['E'][:,i],c=grpF1r['m'][:,i],vmin=-1,vmax=1, marker='_')
panel.set_xlabel(r'$q/k_L$')
panel.set_ylabel(r'$E/E_L$')
plt.colorbar(d)'''
grpF1r.attrs['makeLowestBandFigure']= '''fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)   
pan2.set_title('Lowest band')    
for i in range(3):    
    pan2.scatter(grpF1r['kList'],[i for j in range(grpF1r['kList'].size)],c=grpF1r['pops'][:,i],cmap='Blues',vmin=0.0,vmax=1.0,marker='_',linewidths=10)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'$q/k_L$')'''


datafile3=np.load('SynDimBandStructure_F7_n13_Chern1Flat.npz')
GroupName3='BandStructure_F7'
docstring3=r'''Data for band structure plot as well as populations in the lowest band.
Here, the Raman coupling is flat accross spin states, ie no CG coefficients,
F=7,number of lattice orders included is 27, Raman coupling strength 
$\Omega$=%.3f $E_L$, detuning $\delta$=%.3f, and quadratic shift $\epsilon$=%.3f'''%(datafile3['omega'],datafile3['delta'],datafile3['epsilon'])
grp3 = f.create_group(GroupName3)
grp3.attrs['Note']=docstring3
grp3["kList"] = datafile3['kList']
grp3["kList"].attrs['Note'] = r"List of quasimomenta used in units of $k_L$"
grp3["E"] = datafile3['E'][:,:15]
grp3["E"].attrs['Note'] = r"Energies of the lowest 15 bands in units of $E_L$. Energy in ith band is E[:,i]"
grp3["m"] = datafile3['m'][:,:15]
grp3["m"].attrs['Note'] = r"Magnetizations of the lowest 15 bands. Magnetization in ith band is m[:,i]"
grp3["pops"]= datafile3['pops'][:,0,:]
grp3["pops"].attrs['Note']=r"Lowest band populations in each spin state. Population in each spin state mF is pops[:,mF-7]"
grp3.attrs['makeBandStructureFigure']= '''fig=plt.figure()
panel=fig.add_subplot(1,1,1)
for i in range(15):   
d=panel.scatter(grp3['kList'],grp3['E'][:,i],c=grp3['m'][:,i],vmin=-1,vmax=1, marker='_')
panel.set_xlabel(r'$q/k_L$')
panel.set_ylabel(r'$E/E_L$')
plt.colorbar(d)'''
grp3.attrs['makeLowestBandFigure']= '''fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)   
pan2.set_title('Lowest band')    
for i in range(15):    
    pan2.scatter(grp3['kList'],[i for j in range(grp3['kList'].size)],c=grp3['pops'][:,i],cmap='Blues',vmin=0.0,vmax=1.0,marker='_',linewidths=10)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'$q/k_L$')'''

datafileRfF1=np.load('09Mar2017_Rf_Corrected.npz')
GroupNameRfF1='Data_Rf_F1'
grpRfF1=f.create_group(GroupNameRfF1)
grpRfF1.attrs['Note']=r'''Data for Rf coupled system in F=1, taken on 2017_03_09. 
This is the combined data from kicking in both directions, where the time 
dependence has been corrected by fitting no kick data to a detuning model. Here, Rf coupling strength
was 0.5 $E_L$, lattice depth was 4.4 $E_L$ and quadratic Zeeman shift was 0.02 $E_L$'''
grpRfF1['qlist']=datafileRfF1['qlist']
grpRfF1['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpRfF1['fractionP']=datafileRfF1['fractionP']
grpRfF1['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpRfF1['fraction0']=datafileRfF1['fraction0']
grpRfF1['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpRfF1['fractionM']=datafileRfF1['fractionM']
grpRfF1['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpRfF1['sigmaP']=datafileRfF1['sigmaP']
grpRfF1['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpRfF1['sigma0']=datafileRfF1['sigma0']
grpRfF1['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpRfF1['sigmaM']=datafileRfF1['sigmaM']
grpRfF1['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpRfF1['Mag']=datafileRfF1['fractionP']-datafileRfF1['fractionM']
grpRfF1['Mag'].attrs['Note']=r'Magnetization'
grpRfF1['sigmaMag']=np.sqrt(datafileRfF1['sigmaP']**2.0+datafileRfF1['sigmaM']**2.0)
grpRfF1['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpRfF1.attrs['makePopulationsFigure']='''fig=plt.figure()
pan=fig.add_subplot(1,1,1)
pan.errorbar(grpRfF1['qlist'], grpRfF1['fractionP'],yerr=grpRfF1['sigmaP'],fmt='bo',label='mF=+1')
pan.errorbar(grpRfF1['qlist'], grpRfF1['fraction0'],yerr=grpRfF1['sigma0'],fmt='go',label='mF=0')
pan.errorbar(grpRfF1['qlist'], grpRfF1['fractionM'],yerr=grpRfF1['sigmaM'],fmt='ro',label='mF=-1')
pan.set_xlabel(r'quasimomentum [$k_L$]')
pan.set_ylabel('Corrected fractional populations')'''
grpRfF1.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpRfF1['qlist'],[1 for j in range(grpRfF1['qlist.size'])],c=grpRfF1['fraction'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF1['qlist'],[0 for j in range(grpRfF1['qlist.size'])],c=grpRfF1['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF1['qlist'],[-1 for j in range(grpRfF1['qlist.size'])],c=grpRfF1['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''

datafileRfF2=np.load('22Mar2017_Rf_Corrected.npz')
GroupNameRfF2='Data_Rf_F2'
grpRfF2=f.create_group(GroupNameRfF2)
grpRfF2.attrs['Note']=r'''Data for Rf coupled system in F=2, taken on 2017_03_22. 
This is the combined data from kicking in both directions, where the time 
dependence has been corrected by fitting no kick data to a detuning model. Here, Rf coupling strength
was 0.57 $E_L$, lattice depth was 4.4 $E_L$ and quadratic Zeeman shift was 0.02 $E_L$'''
grpRfF2['qlist']=datafileRfF2['qlist']
grpRfF2['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpRfF2['fractionP2']=datafileRfF2['fractionP2']
grpRfF2['fractionP2'].attrs['Note']=r'fraction in mF=+2 state'
grpRfF2['fractionP']=datafileRfF2['fractionP']
grpRfF2['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpRfF2['fraction0']=datafileRfF2['fraction0']
grpRfF2['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpRfF2['fractionM']=datafileRfF2['fractionM']
grpRfF2['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpRfF2['fractionM2']=datafileRfF2['fractionM2']
grpRfF2['fractionM2'].attrs['Note']=r'fraction in mF=-2 state'
grpRfF2['sigmaP2']=datafileRfF2['sigmaP2']
grpRfF2['sigmaP2'].attrs['Note']=r'Uncertainty in fraction in mF=+2, given by $\sigma/\sqrt(N)$'
grpRfF2['sigmaP']=datafileRfF2['sigmaP']
grpRfF2['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpRfF2['sigma0']=datafileRfF2['sigma0']
grpRfF2['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpRfF2['sigmaM']=datafileRfF2['sigmaM']
grpRfF2['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpRfF2['sigmaM2']=datafileRfF2['sigmaM2']
grpRfF2['sigmaM2'].attrs['Note']=r'Uncertainty in fraction in mF=-2, given by $\sigma/\sqrt(N)$'
grpRfF2['Mag']=2.0*datafileRfF2['fractionP2']+datafileRfF2['fractionP']-datafileRfF2['fractionM']-2.0*datafileRfF2['fractionM2']
grpRfF2['Mag'].attrs['Note']=r'Magnetization'
grpRfF2['sigmaMag']=np.sqrt(4.0*datafileRfF2['sigmaP2']**2.0+datafileRfF2['sigmaP']**2.0+datafileRfF2['sigmaM']**2.0+4.0*datafileRfF2['sigmaM2']**2.0)
grpRfF1['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpRfF2.attrs['makePopulationsFigure']='''fig=plt.figure()
pan=fig.add_subplot(1,1,1)
pan.errorbar(grpRfF2['qlist'], grpRfF2['fractionP'],yerr=grpRfF2['sigmaP'],fmt='bo',label='mF=+1')
pan.errorbar(grpRfF2['qlist'], grpRfF2['fraction0'],yerr=grpRfF2['sigma0'],fmt='go',label='mF=0')
pan.errorbar(grpRfF2['qlist'], grpRfF2['fractionM'],yerr=grpRfF2['sigmaM'],fmt='ro',label='mF=-1')
pan.errorbar(grpRfF2['qlist'], grpRfF2['fractionP2'],yerr=grpRfF2['sigmaP2'],fmt='co',label='mF=+2')
pan.errorbar(grpRfF2['qlist'], grpRfF2['fractionM2'],yerr=grpRfF2['sigmaM2'],fmt='mo',label='mF=-2')
pan.set_xlabel(r'quasimomentum [$k_L$]')
pan.set_ylabel('Corrected fractional populations')'''
grpRfF2.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpRfF2['qlist'],[1 for j in range(grpRfF2['qlist'].size)],c=grpRfF2['fractionP'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF2['qlist'],[0 for j in range(grpRfF2['qlist'].size)],c=grpRfF2['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF2['qlist'],[-1 for j in range(grpRfF2['qlist'].size)],c=grpRfF2['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF2['qlist'],[2 for j in range(grpRfF2['qlist'].size)],c=grpRfF2['fractionP2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpRfF2['qlist'],[-2 for j in range(grpRfF2['qlist'].size)],c=grpRfF2['fractionM2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''

datafileF2_negFlux=np.load('27Feb2017_F2_chern_-1.npz')
GroupNameF2_negFlux='Data_F2_negFlux'
grpF2neg=f.create_group(GroupNameF2_negFlux)
grpF2neg.attrs['Note']=r'''Data and simulation curve for Raman coupled system in F=2
with negative flux, taken on 2017_02_27. 
This is the combined data from kicking in both directions, excluding any data
outside of the first Brillouin zone. Here, Raman coupling strength
was %.2f $E_L$, lattice depth was %.2F $E_L$ and quadratic Zeeman shift was %.2f $E_L$'''%(datafileF2_negFlux['omega'],datafileF2_negFlux['U'],datafileF2_negFlux['epsilon'])
grpF2neg['qlist']=datafileF2_negFlux['qlist']
grpF2neg['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpF2neg['fractionP2']=datafileF2_negFlux['fractionP2']
grpF2neg['fractionP2'].attrs['Note']=r'fraction in mF=+2 state'
grpF2neg['fractionP']=datafileF2_negFlux['fractionP']
grpF2neg['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF2neg['fraction0']=datafileF2_negFlux['fraction0']
grpF2neg['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF2neg['fractionM']=datafileF2_negFlux['fractionM']
grpF2neg['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF2neg['fractionM2']=datafileF2_negFlux['fractionM2']
grpF2neg['fractionM2'].attrs['Note']=r'fraction in mF=-2 state'
grpF2neg['sigmaP2']=datafileF2_negFlux['sigmaP2']
grpF2neg['sigmaP2'].attrs['Note']=r'Uncertainty in fraction in mF=+2, given by $\sigma/\sqrt(N)$'
grpF2neg['sigmaP']=datafileF2_negFlux['sigmaP']
grpF2neg['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpF2neg['sigma0']=datafileF2_negFlux['sigma0']
grpF2neg['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpF2neg['sigmaM']=datafileF2_negFlux['sigmaM']
grpF2neg['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpF2neg['sigmaM2']=datafileF2_negFlux['sigmaM2']
grpF2neg['sigmaM2'].attrs['Note']=r'Uncertainty in fraction in mF=-2, given by $\sigma/\sqrt(N)$'
grpF2neg['qlistSimul']=datafileF2_negFlux['qlistSimul']
grpF2neg['qlistSimul'].attrs['Note']=r'list of quasimomenta for theory curve in units of $k_L$'
grpF2neg['popsSimul']=datafileF2_negFlux['pops']
grpF2neg['popsSimul'].attrs['Note']=r'Fractional populations from theory. Population in each mF state is popsSiml[:,mF+2]'
grpF2neg['Mag']=2.0*datafileF2_negFlux['fractionP2']+datafileF2_negFlux['fractionP']-datafileF2_negFlux['fractionM']-2.0*datafileF2_negFlux['fractionM2']
grpF2neg['Mag'].attrs['Note']=r'Magnetization'
grpF2neg['sigmaMag']=np.sqrt(4.0*datafileF2_negFlux['sigmaP2']**2.0+datafileF2_negFlux['sigmaP']**2.0+datafileF2_negFlux['sigmaM']**2.0+4.0*datafileF2_negFlux['sigmaM2']**2.0)
grpF2neg['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpF2neg.attrs['makePopulationsFigure']='''
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.errorbar(grpF2neg['qlist'],grpF2neg['fractionP'],yerr=grpF2neg['sigmaP'], fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(grpF2neg['qlist'],grpF2neg['fraction0'],yerr=grpF2neg['sigma0'], fmt='go', label=r'$m_F$=0')
pan1.errorbar(grpF2neg['qlist'],grpF2neg['fractionM'],yerr=grpF2neg['sigmaM'], fmt='ro', label=r'$m_F$=-1')
pan1.errorbar(grpF2neg['qlist'],grpF2neg['fractionP2'],yerr=grpF2neg['sigmaP2'], fmt='co', label=r'$m_F$=+2')
pan1.errorbar(grpF2neg['qlist'],grpF2neg['fractionM2'],yerr=grpF2neg['sigmaM2'], fmt='mo', label=r'$m_F$=-2')
for mF in np.arange(-2,3):
    pan1.plot(grpF2neg['qlistSimul'],grpF2neg['popsSimul'][:,mF+S],cDict[mF])
pan1.set_xlabel(r'quasimomentum [$k_L$]')
pan1.set_ylabel('Fractional populations')'''
grpF2neg.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpF2neg['qlist'],[1 for j in range(grpF2neg['qlist'].size)],c=grpF2neg['fractionP'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2neg['qlist'],[0 for j in range(grpF2neg['qlist'].size)],c=grpF2neg['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2neg['qlist'],[-1 for j in range(grpF2neg['qlist'].size)],c=grpF2neg['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2neg['qlist'],[2 for j in range(grpF2neg['qlist'].size)],c=grpF2neg['fractionP2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2neg['qlist'],[-2 for j in range(grpF2neg['qlist'].size)],c=grpF2neg['fractionM2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''


datafileF2_posFlux=np.load('28Feb2017_F2_chern_1.npz')
GroupNameF2_posFlux='Data_F2_posFlux'
grpF2pos=f.create_group(GroupNameF2_posFlux)
grpF2pos.attrs['Note']=r'''Data and simulation curve for Raman coupled system in F=2
with positive flux, taken on 2017_02_28. 
This is the combined data from kicking in both directions, excluding any data
outside of the first Brillouin zone. Here, Raman coupling strength
was %.2f $E_L$, lattice depth was %.2F $E_L$ and quadratic Zeeman shift was %.2f $E_L$'''%(datafileF2_posFlux['omega'],datafileF2_posFlux['U'],datafileF2_posFlux['epsilon'])
grpF2pos['qlist']=datafileF2_posFlux['qlist']
grpF2pos['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpF2pos['fractionP2']=datafileF2_posFlux['fractionP2']
grpF2pos['fractionP2'].attrs['Note']=r'fraction in mF=+2 state'
grpF2pos['fractionP']=datafileF2_posFlux['fractionP']
grpF2pos['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF2pos['fraction0']=datafileF2_posFlux['fraction0']
grpF2pos['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF2pos['fractionM']=datafileF2_posFlux['fractionM']
grpF2pos['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF2pos['fractionM2']=datafileF2_posFlux['fractionM2']
grpF2pos['fractionM2'].attrs['Note']=r'fraction in mF=-2 state'
grpF2pos['sigmaP2']=datafileF2_posFlux['sigmaP2']
grpF2pos['sigmaP2'].attrs['Note']=r'Uncertainty in fraction in mF=+2, given by $\sigma/\sqrt(N)$'
grpF2pos['sigmaP']=datafileF2_posFlux['sigmaP']
grpF2pos['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpF2pos['sigma0']=datafileF2_posFlux['sigma0']
grpF2pos['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpF2pos['sigmaM']=datafileF2_posFlux['sigmaM']
grpF2pos['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpF2pos['sigmaM2']=datafileF2_posFlux['sigmaM2']
grpF2pos['sigmaM2'].attrs['Note']=r'Uncertainty in fraction in mF=-2, given by $\sigma/\sqrt(N)$'
grpF2pos['qlistSimul']=datafileF2_posFlux['qlistSimul']
grpF2pos['qlistSimul'].attrs['Note']=r'list of quasimomenta for theory curve in units of $k_L$'
grpF2pos['popsSimul']=datafileF2_posFlux['pops']
grpF2pos['popsSimul'].attrs['Note']=r'Fractional populations from theory. Population in each mF state is popsSiml[:,mF+2]'
grpF2pos['Mag']=2.0*datafileF2_posFlux['fractionP2']+datafileF2_posFlux['fractionP']-datafileF2_posFlux['fractionM']-2.0*datafileF2_posFlux['fractionM2']
grpF2pos['Mag'].attrs['Note']=r'Magnetization'
grpF2pos['sigmaMag']=np.sqrt(4.0*datafileF2_posFlux['sigmaP2']**2.0+datafileF2_posFlux['sigmaP']**2.0+datafileF2_posFlux['sigmaM']**2.0+4.0*datafileF2_posFlux['sigmaM2']**2.0)
grpF2pos['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpF2pos.attrs['makePopulationsFigure']='''
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.errorbar(grpF2pos['qlist'],grpF2pos['fractionP'],yerr=grpF2pos['sigmaP'], fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(grpF2pos['qlist'],grpF2pos['fraction0'],yerr=grpF2pos['sigma0'], fmt='go', label=r'$m_F$=0')
pan1.errorbar(grpF2pos['qlist'],grpF2pos['fractionM'],yerr=grpF2pos['sigmaM'], fmt='ro', label=r'$m_F$=-1')
pan1.errorbar(grpF2pos['qlist'],grpF2pos['fractionP2'],yerr=grpF2pos['sigmaP2'], fmt='co', label=r'$m_F$=+2')
pan1.errorbar(grpF2pos['qlist'],grpF2pos['fractionM2'],yerr=grpF2pos['sigmaM2'], fmt='mo', label=r'$m_F$=-2')
for mF in np.arange(-2,3):
    pan1.plot(grpF2pos['qlistSimul'],grpF2pos['popsSimul'][:,mF+S],cDict[mF])
pan1.set_xlabel(r'quasimomentum [$k_L$]')
pan1.set_ylabel('Fractional populations')'''
grpF2pos.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpF2pos['qlist'],[1 for j in range(grpF2pos['qlist'].size)],c=grpF2pos['fractionP'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2pos['qlist'],[0 for j in range(grpF2pos['qlist'].size)],c=grpF2pos['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2pos['qlist'],[-1 for j in range(grpF2pos['qlist'].size)],c=grpF2pos['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2pos['qlist'],[2 for j in range(grpF2pos['qlist'].size)],c=grpF2pos['fractionP2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF2pos['qlist'],[-2 for j in range(grpF2pos['qlist'].size)],c=grpF2pos['fractionM2'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''

datafileF1_negFlux=np.load('07Mar2017_F1_chern_-1.npz')
GroupNameF1_negFlux='Data_F1_negFlux'
grpF1neg=f.create_group(GroupNameF1_negFlux)
grpF1neg.attrs['Note']=r'''Data and simulation curve for Raman coupled system in F=1
with negative flux, taken on 2017_03_07. 
This is the combined data from kicking in both directions, excluding any data
outside of the first Brillouin zone. Here, Raman coupling strength
was %.2f $E_L$, lattice depth was %.2F $E_L$ and quadratic Zeeman shift was %.2f $E_L$'''%(datafileF1_negFlux['omega'],datafileF1_negFlux['U'],datafileF1_negFlux['epsilon'])
grpF1neg['qlist']=datafileF1_negFlux['qlist']
grpF1neg['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpF1neg['fractionP']=datafileF1_negFlux['fractionP']
grpF1neg['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF1neg['fraction0']=datafileF1_negFlux['fraction0']
grpF1neg['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF1neg['fractionM']=datafileF1_negFlux['fractionM']
grpF1neg['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF1neg['sigmaP']=datafileF1_negFlux['sigmaP']
grpF1neg['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpF1neg['sigma0']=datafileF1_negFlux['sigma0']
grpF1neg['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpF1neg['sigmaM']=datafileF1_negFlux['sigmaM']
grpF1neg['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpF1neg['qlistSimul']=datafileF1_negFlux['qlistSimul']
grpF1neg['qlistSimul'].attrs['Note']=r'list of quasimomenta for theory curve in units of $k_L$'
grpF1neg['popsSimul']=datafileF1_negFlux['pops']
grpF1neg['popsSimul'].attrs['Note']=r'Fractional populations from theory. Population in each mF state is popsSiml[:,mF+1]'
grpF1neg['Mag']=datafileF1_negFlux['fractionP']-datafileF1_negFlux['fractionM']
grpF1neg['Mag'].attrs['Note']=r'Magnetization'
grpF1neg['sigmaMag']=np.sqrt(datafileF1_negFlux['sigmaP']**2.0+datafileF1_negFlux['sigmaM']**2.0)
grpF1neg['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpF1neg.attrs['makePopulationsFigure']='''
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.errorbar(grpF1neg['qlist'],grpF1neg['fractionP'],yerr=grpF1neg['sigmaP'], fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(grpF1neg['qlist'],grpF1neg['fraction0'],yerr=grpF1neg['sigma0'], fmt='go', label=r'$m_F$=0')
pan1.errorbar(grpF1neg['qlist'],grpF1neg['fractionM'],yerr=grpF1neg['sigmaM'], fmt='ro', label=r'$m_F$=-1')
for mF in np.arange(-1,2):
    pan1.plot(grpF1neg['qlistSimul'],grpF1neg['popsSimul'][:,mF+S],cDict[mF])
pan1.set_xlabel(r'quasimomentum [$k_L$]')
pan1.set_ylabel('Fractional populations')'''
grpF1neg.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpF1neg['qlist'],[1 for j in range(grpF1neg['qlist'].size)],c=grpF1neg['fractionP'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF1neg['qlist'],[0 for j in range(grpF1neg['qlist'].size)],c=grpF1neg['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF1neg['qlist'],[-1 for j in range(grpF1neg['qlist'].size)],c=grpF1neg['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''


datafileF1_posFlux=np.load('08Mar2017_F1_chern_1.npz')
GroupNameF1_posFlux='Data_F1_posFlux'
grpF1pos=f.create_group(GroupNameF1_posFlux)
grpF1pos.attrs['Note']=r'''Data and simulation curve for Raman coupled system in F=1
with positive flux, taken on 2017_03_08. 
This is the combined data from kicking in both directions, excluding any data
outside of the first Brillouin zone. Here, Raman coupling strength
was %.2f $E_L$, lattice depth was %.2F $E_L$ and quadratic Zeeman shift was %.2f $E_L$'''%(datafileF1_posFlux['omega'],datafileF1_posFlux['U'],datafileF1_posFlux['epsilon'])
grpF1pos['qlist']=datafileF1_posFlux['qlist']
grpF1pos['qlist'].attrs['Note']=r'quasimomentum in units of $k_L$'
grpF1pos['fractionP']=datafileF1_posFlux['fractionP']
grpF1pos['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF1pos['fraction0']=datafileF1_posFlux['fraction0']
grpF1pos['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF1pos['fractionM']=datafileF1_posFlux['fractionM']
grpF1pos['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF1pos['sigmaP']=datafileF1_posFlux['sigmaP']
grpF1pos['sigmaP'].attrs['Note']=r'Uncertainty in fraction in mF=+1, given by $\sigma/\sqrt(N)$'
grpF1pos['sigma0']=datafileF1_posFlux['sigma0']
grpF1pos['sigma0'].attrs['Note']=r'Uncertainty in fraction in mF=0, given by $\sigma/\sqrt(N)$'
grpF1pos['sigmaM']=datafileF1_posFlux['sigmaM']
grpF1pos['sigmaM'].attrs['Note']=r'Uncertainty in fraction in mF=-1, given by $\sigma/\sqrt(N)$'
grpF1pos['qlistSimul']=datafileF1_posFlux['qlistSimul']
grpF1pos['qlistSimul'].attrs['Note']=r'list of quasimomenta for theory curve in units of $k_L$'
grpF1pos['popsSimul']=datafileF1_posFlux['pops']
grpF1pos['popsSimul'].attrs['Note']=r'Fractional populations from theory. Population in each mF state is popsSiml[:,mF+1]'
grpF1pos['Mag']=datafileF1_posFlux['fractionP']-datafileF1_posFlux['fractionM']
grpF1pos['Mag'].attrs['Note']=r'Magnetization'
grpF1pos['sigmaMag']=np.sqrt(datafileF1_posFlux['sigmaP']**2.0+datafileF1_posFlux['sigmaM']**2.0)
grpF1pos['sigmaMag'].attrs['Note']=r'Propagated uncertainty in magnetization'
grpF1pos.attrs['makePopulationsFigure']='''
cDict={}
cDict[-2]='m-'
cDict[-1]='r-'
cDict[0]='g-'
cDict[1]='b-'
cDict[2]='c-'
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
pan1.errorbar(grpF1pos['qlist'],grpF1pos['fractionP'],yerr=grpF1pos['sigmaP'], fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(grpF1pos['qlist'],grpF1pos['fraction0'],yerr=grpF1pos['sigma0'], fmt='go', label=r'$m_F$=0')
pan1.errorbar(grpF1pos['qlist'],grpF1pos['fractionM'],yerr=grpF1pos['sigmaM'], fmt='ro', label=r'$m_F$=-1')
for mF in np.arange(-1,2):
    pan1.plot(grpF1pos['qlistSimul'],grpF1pos['popsSimul'][:,mF+S],cDict[mF])
pan1.set_xlabel(r'quasimomentum [$k_L$]')
pan1.set_ylabel('Fractional populations')'''
grpF1pos.attrs['makeBarcodePlot']='''
fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.scatter(grpF1pos['qlist'],[1 for j in range(grpF1pos['qlist'].size)],c=grpF1pos['fractionP'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF1pos['qlist'],[0 for j in range(grpF1pos['qlist'].size)],c=grpF1pos['fraction0'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.scatter(grpF1pos['qlist'],[-1 for j in range(grpF1pos['qlist'].size)],c=grpF1pos['fractionM'],vmin=0.0,vmax=1.0,cmap='Blues', marker='_',linewidths=30)
pan2.set_ylabel('Synthetic lattice site')
pan2.set_xlabel(r'Quasimomentum [$k_L$]')'''

datafileImages=np.load('dataForImagesFigure.npz')
GroupNameImages='F2images'
grpF2images=f.create_group(GroupNameImages)
grpF2images.attrs['Note']=r'''Sample data for kicking both left and right in F=2, taken on 27th of February 2017. 
There are two arrays of images - one raw, divided into separate spin states and one numerically 
'demapped'. Both arrays have the first index as the image number, listed in order of quasimomentum, 
so the first 15 are from a negative kick data (starting from the longest kick time) and the last 
15 are from the positive kick data (starting from shortest time). '''
grpF2images['odF2']=datafileImages['odF2']
grpF2images['odF2'].attrs['Note']='''array of raw images for both kick directions split up
into different spin states. The indexing is as follows: odF2[file,spin,ypixels,xpixels]. File goes from 0 to 29,
with 0-14 being negative kick data and 15-29 being positive kick data. Spin goes from 0 to 4, with 0 representing 
mF=-2 and 4 representing mF=2. ypixels goes from 0 to 309 (y direction in images, along the lattice direction) 
and xpixels goes from 0 to 59 (x direction in images, perpendicular to lattice direction)'''
grpF2images['odF2demap']=datafileImages['odF2demap']
grpF2images['odF2demap'].attrs['Note']='''array of 'demapped' images for both kick directions split up
into different spin states. The indexing is as follows: odF2[file,spin,quasimomentum,xpixels]. File goes from 0 to 29,
with 0-14 being negative kick data and 15-29 being positive kick data. Spin goes from 0 to 4, with 0 representing 
mF=-2 and 4 representing mF=2. quasimomentum goes from 0 to 61, with zero representing -k_L and 61 representing +k_L,
and xpixels goes from 0 to 59 (x direction in images, perpendicular to lattice direction)'''
grpF2images['qlist']=datafileImages['qlist']
grpF2images['qlist'].attrs['Note']='''Quasimomentum (as measured by taking the location of the central order) of each
image file, in the same order as the first index of the image arrays - ie first one is the longest hold time from the 
negative kick data (most negative q) and last one is the longest hold time of the positive kick data (most positive q)'''
grpF2images['tlist']=datafileImages['tlist']
grpF2images['tlist'].attrs['Note']='''Actual kick times for each image, in the same order as the quasimomentum
and the first index in the image arrays'''

datafilePulsing1=np.load('22Aug2017_files_4-33.npz')
GroupNamePulsing1='PulsingData1'
grpF2pulsing1=f.create_group(GroupNamePulsing1)
grpF2pulsing1.attrs['Note']=r'''Sample pulsing data from 22Aug2017. The lattice depth here was 4.4 E_L and the 
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

datafilePulsing2=np.load('22Aug2017_files_34-63.npz')
GroupNamePulsing2='PulsingData2'
grpF2pulsing2=f.create_group(GroupNamePulsing2)
grpF2pulsing2.attrs['Note']=r'''Sample pulsing data from 22Aug2017. The lattice depth here was 4.4 E_L and the 
fitted Raman coupling strenth was 0.52 E_L with a fitted detuning of 0.002 E_L. Quadratic zeeman shift was 0.02 E_L
The lattice was turned on adiabatically, and the RAman was turned on in 300 us'''
grpF2pulsing2['tlist']=datafilePulsing2['tlist']
grpF2pulsing2['tlist'].attrs['Note']=r'Raman pulse time in seconds'
grpF2pulsing2['fractionP2']=datafilePulsing2['fractionP2']
grpF2pulsing2['fractionP2'].attrs['Note']=r'fraction in mF=+2 state'
grpF2pulsing2['fractionP']=datafilePulsing2['fractionP']
grpF2pulsing2['fractionP'].attrs['Note']=r'fraction in mF=+1 state'
grpF2pulsing2['fraction0']=datafilePulsing2['fraction0']
grpF2pulsing2['fraction0'].attrs['Note']=r'fraction in mF=0 state'
grpF2pulsing2['fractionM']=datafilePulsing2['fractionM']
grpF2pulsing2['fractionM'].attrs['Note']=r'fraction in mF=-1 state'
grpF2pulsing2['fractionM2']=datafilePulsing2['fractionM2']
grpF2pulsing2['fractionM2'].attrs['Note']=r'fraction in mF=-2 state'

print('Size of %s is %.2fMiB' % (h5name, float(os.stat(h5name).st_size)/1024**2)) 


f.close()