# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:23:27 2017

@author: dng5
"""

import numpy as np

import matplotlib.pyplot as plt
from scipy import optimize
import lineFit

#dataFileList =np.array(['19Oct2017_files_34-303_dummy_10.0.npz',
#                        '19Oct2017_files_34-303_dummy_3.0.npz',
#                        '19Oct2017_files_34-303_dummy_8.0.npz',
#                        '19Oct2017_files_34-303_dummy_2.0.npz',
#                        '19Oct2017_files_34-303_dummy_7.0.npz',
#                        '19Oct2017_files_34-303_dummy_9.0.npz',
#                        '19Oct2017_files_34-303_dummy_5.0.npz',
#                        '19Oct2017_files_34-303_dummy_4.0.npz',
#                        '19Oct2017_files_34-303_dummy_6.0.npz'])
#saveFileName='12Apr2017_F2Calib_newBgnd'
#dataFileList =np.array(['20Oct2017_files_4-153_dummy_1.0.npz',
#                        '20Oct2017_files_4-153_dummy_2.0.npz',
#                        '20Oct2017_files_4-153_dummy_3.0.npz',
#                        '20Oct2017_files_4-153_dummy_4.0.npz',
#                        '20Oct2017_files_4-153_dummy_5.0.npz',
#                        '20Oct2017_files_4-153_dummy_6.0.npz',
#                        '20Oct2017_files_4-153_dummy_7.0.npz',
#                        '20Oct2017_files_4-153_dummy_8.0.npz',
#                        '20Oct2017_files_4-153_dummy_9.0.npz',
#                        '20Oct2017_files_4-153_dummy_10.0.npz'])
dataFileList=np.array(['02Nov2017_files_4-113_dummy_1.0.npz',
                       '02Nov2017_files_4-113_dummy_2.0.npz',
                       '02Nov2017_files_4-113_dummy_3.0.npz',
                       '02Nov2017_files_4-113_dummy_4.0.npz',
                       '02Nov2017_files_4-113_dummy_5.0.npz',
                       '02Nov2017_files_4-113_dummy_6.0.npz',
                       '02Nov2017_files_4-113_dummy_7.0.npz',
                       '02Nov2017_files_4-113_dummy_8.0.npz',
                       '02Nov2017_files_4-113_dummy_9.0.npz',
                       '02Nov2017_files_4-113_dummy_10.0.npz'])

tAll=[]

fP=[]
fM=[]
f0=[] 


numP=[]
num0=[]
numM=[]


q=[]
odFiltAll=[]
fileNums=[]

imbal=[]

cutoff = 1.0
for fileName in dataFileList:
    
    dataFile=np.load(fileName)

    times=dataFile['tlist']
    sort=np.argsort(times)

    
    numPloc=dataFile['numP'][sort]
    num0loc=dataFile['num0'][sort]
    numMloc=dataFile['numM'][sort]
        
    kImbal = np.abs(numPloc-numMloc)/(numPloc+numMloc)
    
    imbalGood = (kImbal<cutoff)
    
    times=times[sort][imbalGood]
    qlist=dataFile['qlist']
    qlist=qlist[sort][imbalGood]
                
    fractionP=dataFile['fractionP'][sort][imbalGood]
    fraction0=dataFile['fraction0'][sort][imbalGood]
    fractionM=dataFile['fractionM'][sort][imbalGood]
    
    numPloc=dataFile['numP'][sort][imbalGood]
    num0loc=dataFile['num0'][sort][imbalGood]
    numMloc=dataFile['numM'][sort][imbalGood]
    kImbal=kImbal[imbalGood]

    for i,time in enumerate(times):
        try:
            ind = tAll.index(time)

            fP[ind].append(fractionP[i])
            f0[ind].append(fraction0[i])
            fM[ind].append(fractionM[i])
            q[ind].append(qlist[i])

            numP[ind].append(numPloc[i])
            num0[ind].append(num0loc[i])
            numM[ind].append(numMloc[i])
            imbal[ind].append(kImbal[i])

        except ValueError:
            tAll.append(time)

            fP.append([fractionP[i]])
            f0.append([fraction0[i]])
            fM.append([fractionM[i]])

            numP.append([numPloc[i]])
            num0.append([num0loc[i]])
            numM.append([numMloc[i]])

               
            q.append([qlist[i]])
            imbal.append([kImbal[i]])


   
fractionP=np.array([np.average(fP[i]) for i in range(len(fP))]) 
sigmaP= np.array([np.std(fP[i]) for i in range(len(fP))]) 
fraction0=np.array([np.average(f0[i]) for i in range(len(f0))]) 
sigma0=np.array([np.std(f0[i]) for i in range(len(f0))]) 
fractionM=np.array([np.average(fM[i]) for i in range(len(fM))])  
sigmaM=np.array([np.std(fM[i]) for i in range(len(fM))])


qlist=np.array([np.average(q[i]) for i in range(len(q))])  
sigmaQ=np.array([np.std(q[i]) for i in range(len(q))])
time=np.array(tAll)   
numPoints=np.array([len(q[i]) for i in range(len(q))])



fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(time,fractionP,yerr=sigmaP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(time,fraction0,yerr=sigma0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
pan1.errorbar(time,fractionM,yerr=sigmaM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
#pan1.plot(time,gauss(time,0.00006,0.15,0.000008,0.0),'b.')

pan1.set_xlabel('Time [s]')
pan1.set_ylabel('Fractional populations')

#x0,A,sigma,offset = gaussFit(time,fraction0,'pulse time [s]','fraction in k = 0',p0=(0.00006,-0.15,0.00001,1.0))


fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(time,fractionP,yerr=sigmaP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(time,fraction0,yerr=sigma0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
pan1.errorbar(time,fractionM,yerr=sigmaM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
#tforfit=np.linspace(np.min(time),np.max(time),100)
#pan1.plot(tforfit,gauss(tforfit,x0,A,sigma,offset),'g-')

pan1.set_xlabel('Time [s]')
pan1.set_ylabel('Fractional populations')
#np.savez(saveFileName,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,
#         qlist=qlist,sigmaQ=sigmaQ,tlist=time,sigmaP=sigmaP,sigma0=sigma0,sigmaM=sigmaM,numPoints=numPoints,
#         mag=mag,sigmaMag=sigmaMag,qlistFit=qlistFit,newQlist=newQlist,newP=newP,
#         newSigmaP=newSigmaP,new0=new0,newSigma0=newSigma0,newM=newM,newSigmaM=newSigmaM)

A,f,phi,offset,resids,chi2,uncerts=sineFit(time,qlist,'hold time','measured q',p0=(0.06,28,0,-0.1))
