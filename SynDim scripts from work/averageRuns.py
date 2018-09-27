# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 17:20:49 2016

@author: dng5
"""
import numpy as np

import matplotlib.pyplot as plt
from scipy import optimize
import lineFit

getNums=False

def bestFitNewX(newXlist,xlist,ylist,sigmaXlist,sigmaYlist):
    nNew=newXlist.size
    nOld=xlist.size
    xDiffMat= (np.array([xlist]*nNew)-np.array([newXlist]*nOld).transpose())**2.0
    yMat=np.array([ylist]*nNew)
    sigmaXmat=np.array([sigmaXlist]*nNew)
    sigmaYmat=np.array([sigmaYlist]*nNew)
    A=np.sum(yMat*np.exp(-xDiffMat/(2.0*sigmaXmat**2.0))/(sigmaYmat**2.0),axis=1)
    B=np.sum(np.exp(-xDiffMat/(2.0*sigmaXmat**2.0))/(sigmaYmat**2.0),axis=1)
    C=np.sum(np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    D=np.sum(yMat*np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    E=np.sum((yMat**2.0)*np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    F=np.sum((yMat**2.0)*(xDiffMat**2.0/(sigmaXmat**4.0))*np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    G=np.sum(yMat*(xDiffMat**2.0/(sigmaXmat**4.0))*np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    H=np.sum((xDiffMat**2.0/(sigmaXmat**4.0))*np.exp(-xDiffMat/(sigmaXmat**2.0))/(sigmaYmat**4.0),axis=1)
    newYlist=A/B
    newSigmaYlist=np.sqrt((F+4.0*E)*B**2.0-2.0*A*B*(G+4.0*D)+(H+4.0*C)*A**2.0)/(B**2.0)
    
    return newYlist,newSigmaYlist
#dataFileList =np.array(['13Dec2016_files_254-303.npz','13Dec2016_files_304-353.npz','13Dec2016_files_454-503_v2.npz',
#                        '13Dec2016_files_554-603_v2.npz','13Dec2016_files_654-703_v2.npz'])
#
#saveFileName='13Dec2016_posKick'

#
#dataFileList =np.array(['30Jan2017_files_184-213.npz','30Jan2017_files_274-303.npz','30Jan2017_files_364-393.npz',
#                        '31Jan2017_files_94-123.npz','31Jan2017_files_184-213.npz','31Jan2017_files_274-303.npz',
#                        '01Feb2017_files_94-122.npz','01Feb2017_files_184-213.npz','01Feb2017_files_274-303.npz'])
#saveFileName='rfF2_negKick'  
#S=2                      
#checkCoh=True
#kick='neg'
#dataFileList =np.array(['30Jan2017_files_154-183.npz','30Jan2017_files_244-273.npz','30Jan2017_files_334-363.npz',
#                        '31Jan2017_files_64-93.npz','31Jan2017_files_154-183.npz','31Jan2017_files_244-273.npz',
#                        '01Feb2017_files_64-93.npz','01Feb2017_files_154-183.npz','01Feb2017_files_244-273.npz'])
#saveFileName='rfF2_posKick'
#S=2                      
#checkCoh=True
#kick='pos'
#dataFileList =np.array(['30Jan2017_files_124-153.npz','30Jan2017_files_215-243.npz','30Jan2017_files_304-333.npz',
#                        '31Jan2017_files_34-63.npz','31Jan2017_files_124-153.npz','31Jan2017_files_214-243.npz',
#                        '01Feb2017_files_34-63.npz','01Feb2017_files_124-153.npz','01Feb2017_files_214-243.npz'])
#saveFileName='rfF2_noKick'
#S=2                      
#checkCoh=True
#kick='neg'
#dataFileList =np.array(['27Feb2017_files_124-153.npz','27Feb2017_files_184-213.npz','27Feb2017_files_244-273.npz',
#                        '27Feb2017_files_304-333.npz','27Feb2017_files_364-393.npz','27Feb2017_files_424-453.npz'])
#saveFileName='27Feb2017_posKick'
#S=2
#kick='pos'
#checkCoh=False
#dataFileList =np.array(['27Feb2017_files_154-183.npz','27Feb2017_files_214-243.npz','27Feb2017_files_274-303.npz',
#                        '27Feb2017_files_334-363.npz','27Feb2017_files_394-423.npz','27Feb2017_files_454-483.npz'])
#saveFileName='27Feb2017_negKick'
#S=2
#kick='neg'
#checkCoh=False
#dataFileList =np.array(['28Feb2017_files_394-423.npz','28Feb2017_files_214-243.npz','28Feb2017_files_274-303.npz',
#                        '28Feb2017_files_334-363.npz','28Feb2017_files_454-483.npz','28Feb2017_files_514-543.npz'])
#saveFileName='28Feb2017_posKick'
#S=2
#kick='posPlus'
#checkCoh=False

#dataFileList =np.array(['28Feb2017_files_544-573.npz','28Feb2017_files_244-273.npz','28Feb2017_files_304-333.npz',
#                        '28Feb2017_files_364-393.npz','28Feb2017_files_424-453.npz','28Feb2017_files_484-513.npz'])
#saveFileName='28Feb2017_negKick'
#S=2
#kick='negPlus'
#checkCoh=False

#dataFileList =np.array(['01Mar2017_files_394-423.npz','01Mar2017_files_424-453.npz','01Mar2017_files_454-483.npz',
#                        '01Mar2017_files_484-513.npz','01Mar2017_files_514-543.npz','01Mar2017_files_544-573.npz',])
#saveFileName='01Mar2017_secondBand_negKick'

dataFileList =np.array(['07Mar2017_files_304-333.npz','07Mar2017_files_364-393.npz','07Mar2017_files_424-453.npz',
                        '07Mar2017_files_484-513.npz','07Mar2017_files_544-573.npz','07Mar2017_files_604-633.npz'])
saveFileName='07Mar2017_negKick'
S=1
kick='neg'
checkCoh=False

#dataFileList =np.array(['07Mar2017_files_274-303.npz','07Mar2017_files_334-363.npz','07Mar2017_files_454-483.npz',
#                        '07Mar2017_files_514-543.npz','07Mar2017_files_574-603.npz','07Mar2017_files_634-663.npz'])
#saveFileName='07Mar2017_posKick'
#S=1
#kick='newrule'
#checkCoh=False

#dataFileList =np.array(['07Mar2017_files_724-753.npz','08Mar2017_files_64-93.npz','08Mar2017_files_124-153.npz',
#                        '08Mar2017_files_184-213.npz','08Mar2017_files_244-273.npz','08Mar2017_files_304-333.npz'])
#saveFileName='08Mar2017_posKick'
#kick='pos'
#S=1
#checkCoh=False

#dataFileList =np.array(['07Mar2017_files_754-783.npz','08Mar2017_files_94-123.npz','08Mar2017_files_154-183.npz',
#                        '08Mar2017_files_214-243.npz','08Mar2017_files_274-303.npz','08Mar2017_files_334-363.npz'])
#saveFileName='08Mar2017_negKick'
#checkCoh=False
#S=1
#kick='neg'
##
#dataFileList =np.array(['09Mar2017_files_274-303.npz','09Mar2017_files_364-393.npz','09Mar2017_files_454-483.npz',
#                        '09Mar2017_files_544-573.npz','09Mar2017_files_634-663.npz','09Mar2017_files_724-753.npz'])
#saveFileName='09Mar2017_noKick'
#S=1
#kick='no'
#checkCoh=True
##
#dataFileList =np.array(['09Mar2017_files_304-333.npz','09Mar2017_files_394-423.npz','09Mar2017_files_484-513.npz',
#                        '09Mar2017_files_574-603.npz','09Mar2017_files_664-693.npz','09Mar2017_files_754-783.npz'])
#saveFileName='09Mar2017_posKick'
#S=1
#kick='pos'
#checkCoh=True
#
#dataFileList =np.array(['09Mar2017_files_244-273.npz','09Mar2017_files_334-363.npz','09Mar2017_files_424-453.npz',
#                        '09Mar2017_files_514-543.npz','09Mar2017_files_604-633.npz','09Mar2017_files_694-723.npz'])
#saveFileName='09Mar2017_negKick'
#S=1
#kick='neg'
#checkCoh=True
#
#dataFileList =np.array(['22Mar2017_files_34-63.npz','22Mar2017_files_124-153.npz','22Mar2017_files_214-243.npz',
#                        '22Mar2017_files_304-330.npz','22Mar2017_files_394-423.npz','22Mar2017_files_484-513.npz'])
#saveFileName='22Mar2017_negKick'
#S=2
#kick='neg'
#checkCoh=True
#
#dataFileList =np.array(['22Mar2017_files_64-93.npz','22Mar2017_files_154-183.npz','22Mar2017_files_244-273.npz',
#                        '22Mar2017_files_334-363.npz','22Mar2017_files_424-453.npz','22Mar2017_files_514-543.npz'])
#saveFileName='22Mar2017_posKick'
#S=2
#kick='pos'
#checkCoh=True

#dataFileList =np.array(['22Mar2017_files_94-123.npz','22Mar2017_files_184-213.npz','22Mar2017_files_274-303.npz',
#                        '22Mar2017_files_364-393.npz','22Mar2017_files_454-483.npz','22Mar2017_files_544-573.npz'])
#saveFileName='22Mar2017_noKick'
#S=2
#kick='no'
#checkCoh=True

#dataFileList =np.array(['10Apr2017_files_4-33.npz','10Apr2017_files_64-93.npz','10Apr2017_files_94-123.npz',
#                        '10Apr2017_files_124-153.npz','10Apr2017_files_154-183.npz','10Apr2017_files_184-213.npz','10Apr2017_files_214-218.npz'])
#saveFileName='10Apr2017_F1Calib_newBgnd'
#S=1
#kick='no'
#checkCoh=False
#getNums=True

#dataFileList =np.array(['10Apr2017_files_249-278.npz','10Apr2017_files_279-308.npz','10Apr2017_files_309-338.npz',
#                        '10Apr2017_files_339-368.npz','10Apr2017_files_370-399.npz','10Apr2017_files_400-429.npz'])
#saveFileName='10Apr2017_F2Calib_newBgnd'
#S=2
#kick='no'
#checkCoh=False
#getNums=True


#dataFileList =np.array(['12Apr2017_files_4-33.npz','12Apr2017_files_34-63.npz','12Apr2017_files_64-93.npz',
#                        '12Apr2017_files_94-123.npz','12Apr2017_files_124-153.npz','12Apr2017_files_154-183.npz'])
#saveFileName='12Apr2017_F2Calib_newBgnd'
#S=2
#kick='no'
#checkCoh=False
#getNums=True


#dataFileList =np.array(['12Apr2017_files_184-213.npz'])
#saveFileName='12Apr2017_F2Calib2'
#S=2
#kick='no'
#checkCoh=False
#getNums=True
tAll=[]
fP2=[]
fP=[]
fM=[]
f0=[] 
fM2=[] 
numP2=[]
numP=[]
num0=[]
numM=[]
numM2=[]
num1=[]
num2=[]
numP2norm=[]
numPnorm=[]
num0norm=[]
numMnorm=[]
numM2norm=[]
q=[]
odFiltAll=[]
fileNums=[]


for fileName in dataFileList:
    
    dataFile=np.load(fileName)
    imbal=dataFile['imbalArray']
    signalGood=dataFile['signalGood']
    cutoff=0.2
    fieldGoodArray=((np.abs(imbal)<cutoff) & (signalGood))
    fileGood=(fieldGoodArray)
    if checkCoh:
        coherent=dataFile['coherentList']
        fileGood=(fieldGoodArray & coherent)
    times=dataFile['tlist'][fileGood]
    sort=np.argsort(times)
    times=times[sort]
    qlist=dataFile['qlist'][fileGood]
    qlist=qlist[sort]
    if kick=='neg':
        for i in range(qlist.size-1):
            if qlist[i+1]<qlist[i]-1.5:
                qlist[i+1]=qlist[i+1]+2.0
#                
    if kick=='pos':
        for i in range(qlist.size-1):
            if qlist[i+1]>qlist[i]+1.5:
                qlist[i+1]=qlist[i+1]-2.0
                
    if kick=='newrule':
        for i in range(qlist.size):
            if qlist[i]>-1.65:
                qlist[i]=qlist[i]-2.0
    if kick=='negPlus':
        if qlist[0]>-1.5:
            qlist[0]=qlist[0]-2.0
        for i in range(qlist.size-1):
            if qlist[i+1]<qlist[i]-1.5:
                qlist[i+1]=qlist[i+1]+2.0
                
    if kick=='posPlus':
        if qlist[0]<-1.0:
            qlist[0]=qlist[0]+2.0
        for i in range(qlist.size-1):
            if qlist[i+1]>qlist[i]+1.5:
                qlist[i+1]=qlist[i+1]-2.0
            elif qlist[i+1]<qlist[i]-1.5:
                qlist[i+1]=qlist[i+1]+2
                
                
    fractionP=dataFile['fractionP'][fileGood][sort]
    fraction0=dataFile['fraction0'][fileGood][sort]
    fractionM=dataFile['fractionM'][fileGood][sort]
    
    if getNums:
        num1loc=dataFile['num1Array'][fileGood][sort]
        num2loc=dataFile['num2Array'][fileGood][sort]
        
        numPloc=dataFile['numP'][fileGood][sort]
        num0loc=dataFile['num0'][fileGood][sort]
        numMloc=dataFile['numM'][fileGood][sort]
        
        numPlocN=numPloc*2.0/(num1loc+num2loc)
        num0locN=num0loc*2.0/(num1loc+num2loc)
        numMlocN=numMloc*2.0/(num1loc+num2loc)
        
    if S==2:
        fractionP2=dataFile['fractionP2'][fileGood][sort]
        fractionM2=dataFile['fractionM2'][fileGood][sort]
        if getNums:
            numP2loc=dataFile['numP2'][fileGood][sort]
            numM2loc=dataFile['numM2'][fileGood][sort]
            numP2locN=numP2loc*2.0/(num1loc+num2loc)
            numM2locN=numM2loc*2.0/(num1loc+num2loc)
        
        

#    num1loc=dataFile['num1'][fieldGoodArray][sort]
#    num2loc=dataFile['num2'][fieldGoodArray][sort]
#    num3loc=dataFile['num3'][fieldGoodArray][sort]

    for i,time in enumerate(times):
        try:
            ind = tAll.index(time)
            if S==2:
                fP2[ind].append(fractionP2[i])
                fM2[ind].append(fractionM2[i])
                if getNums:
                    numP2[ind].append(numP2loc[i])
                    numM2[ind].append(numM2loc[i])
                    numP2norm[ind].append(numP2locN[i])
                    numM2norm[ind].append(numM2locN[i])
            fP[ind].append(fractionP[i])
            f0[ind].append(fraction0[i])
            fM[ind].append(fractionM[i])
            q[ind].append(qlist[i])
            if getNums:
                numP[ind].append(numPloc[i])
                num0[ind].append(num0loc[i])
                numM[ind].append(numMloc[i])
                numPnorm[ind].append(numPlocN[i])
                num0norm[ind].append(num0locN[i])
                numMnorm[ind].append(numMlocN[i])
                num1[ind].append(num1loc[i])
                num2[ind].append(num2loc[i])
        except ValueError:
            tAll.append(time)
            if S==2:
                fP2.append([fractionP2[i]])
                fM2.append([fractionM2[i]])
                if getNums:
                    numP2.append([numP2loc[i]])
                    numM2.append([numM2loc[i]])
                    numP2norm.append([numP2locN[i]])
                    numM2norm.append([numM2locN[i]])
            fP.append([fractionP[i]])
            f0.append([fraction0[i]])
            fM.append([fractionM[i]])
            
            if getNums:
                numP.append([numPloc[i]])
                num0.append([num0loc[i]])
                numM.append([numMloc[i]])
                numPnorm.append([numPlocN[i]])
                num0norm.append([num0locN[i]])
                numMnorm.append([numMlocN[i]])
                num1.append([num1loc[i]])
                num2.append([num2loc[i]])
               
            q.append([qlist[i]])
#            num1.append([num1loc[i]])
#            num2.append([num2loc[i]])
#            num3.append([num3loc[i]])
            
            
if S==2:
    fractionP2=np.array([np.average(fP2[i]) for i in range(len(fP2))]) 
    sigmaP2= np.array([np.std(fP2[i]) for i in range(len(fP2))]) 
    fractionM2=np.array([np.average(fM2[i]) for i in range(len(fM2))])  
    sigmaM2=np.array([np.std(fM2[i]) for i in range(len(fM2))])
    if getNums:
        nM2=np.array([np.average(numM2[i]) for i in range(len(numM2))])  
        signM2=np.array([np.std(numM2[i]) for i in range(len(numM2))])
        nP2=np.array([np.average(numP2[i]) for i in range(len(numP2))])  
        signP2=np.array([np.std(numP2[i]) for i in range(len(numP2))])
        nM2norm=np.array([np.average(numM2norm[i]) for i in range(len(numM2norm))])  
        signM2norm=np.array([np.std(numM2norm[i]) for i in range(len(numM2norm))])
        nP2norm=np.array([np.average(numP2norm[i]) for i in range(len(numP2norm))])  
        signP2norm=np.array([np.std(numP2norm[i]) for i in range(len(numP2norm))])
   
fractionP=np.array([np.average(fP[i]) for i in range(len(fP))]) 
sigmaP= np.array([np.std(fP[i]) for i in range(len(fP))]) 
fraction0=np.array([np.average(f0[i]) for i in range(len(f0))]) 
sigma0=np.array([np.std(f0[i]) for i in range(len(f0))]) 
fractionM=np.array([np.average(fM[i]) for i in range(len(fM))])  
sigmaM=np.array([np.std(fM[i]) for i in range(len(fM))])
if getNums:
    nM=np.array([np.average(numM[i]) for i in range(len(numM))])  
    signM=np.array([np.std(numM[i]) for i in range(len(numM))])
    nP=np.array([np.average(numP[i]) for i in range(len(numP))])  
    signP=np.array([np.std(numP[i]) for i in range(len(numP))])
    n0=np.array([np.average(num0[i]) for i in range(len(num0))])  
    sign0=np.array([np.std(num0[i]) for i in range(len(num0))])
    nMnorm=np.array([np.average(numMnorm[i]) for i in range(len(numMnorm))])  
    signMnorm=np.array([np.std(numMnorm[i]) for i in range(len(numMnorm))])
    nPnorm=np.array([np.average(numPnorm[i]) for i in range(len(numPnorm))])  
    signPnorm=np.array([np.std(numPnorm[i]) for i in range(len(numPnorm))])
    n0norm=np.array([np.average(num0norm[i]) for i in range(len(num0norm))])  
    sign0norm=np.array([np.std(num0norm[i]) for i in range(len(num0norm))])
    n1=np.array([np.average(num1[i]) for i in range(len(num1))])  
    sigma1=np.array([np.std(num1[i]) for i in range(len(num1))])
    n2=np.array([np.average(num2[i]) for i in range(len(num2))])  
    sigma2=np.array([np.std(num2[i]) for i in range(len(num2))])

qlist=np.array([np.average(q[i]) for i in range(len(q))])  
sigmaQ=np.array([np.std(q[i]) for i in range(len(q))])
time=np.array(tAll)   
numPoints=np.array([len(q[i]) for i in range(len(q))])
#num1Tot=np.array([np.average(num1[i]) for i in range(len(num1))])  
#sigma1Tot=np.array([np.std(num1[i]) for i in range(len(num1))])
#num2Tot=np.array([np.average(num2[i]) for i in range(len(num2))])  
#sigma2Tot=np.array([np.std(num2[i]) for i in range(len(num2))])
#num3Tot=np.array([np.average(num3[i]) for i in range(len(num3))])  
#sigma3Tot=np.array([np.std(num3[i]) for i in range(len(num3))])
#np.savez(saveFileName,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,
#         qlist=qlist,sigmaQ=sigmaQ,tlist=time,sigmaP=sigmaP,sigma0=sigma0,sigmaM=sigmaM,numPoints=numPoints,
#         num1Tot=num1Tot,num2Tot=num2Tot,num3Tot=num3Tot, sigma1Tot=sigma1Tot,
#         sigma2Tot=sigma2Tot,sigma3Tot=sigma3Tot)
mag=fractionP-fractionM#-2.0*fractionM2+2.0*fractionP2
sigmaMag=np.sqrt(sigmaP**2.0+sigmaM**2.0)/np.sqrt(numPoints)
if S==2: 
    mag=fractionP-fractionM-2.0*fractionM2+2.0*fractionP2
    sigmaMag=np.sqrt(4.0*sigmaP2+sigmaP**2.0+sigmaM**2.0+4.0*sigmaM2)/np.sqrt(numPoints)

#(A,x0,Gamma), cov=optimize.curve_fit(lor,time,fraction0)
#xrangefit=np.linspace(np.min(time),np.max(time),600)
#data_fitted=lor(xrangefit,*(A,x0,Gamma))
    

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(time,fractionP,yerr=sigmaP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(time,fraction0,yerr=sigma0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
pan1.errorbar(time,fractionM,yerr=sigmaM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan1.errorbar(time,fractionP2,yerr=sigmaP2/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
    pan1.errorbar(time,fractionM2,yerr=sigmaM2/np.sqrt(numPoints), fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan1.set_xlabel('Time [s]')
pan1.set_ylabel('Fractional populations')
legend()

if getNums:
    fig2=plt.figure()
    pan2=fig2.add_subplot(1,1,1)
    pan2.errorbar(time,nP,yerr=signP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
    pan2.errorbar(time,n0,yerr=sign0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
    pan2.errorbar(time,nM,yerr=signM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
    tot=nM+nP+n0
    sigmaTot=np.sqrt(signM**2.0+sign0**2.0+signP**2.0)
    if S==2:
        pan2.errorbar(time,nP2,yerr=signP2/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
        pan2.errorbar(time,nM2,yerr=signM2/np.sqrt(numPoints), fmt='mo', label=r'$m_F$=-2')
        tot=nM+nP+n0+nP2+nM2
        sigmaTot=np.sqrt(signM**2.0+sign0**2.0+signP**2.0+signP2**2.0+signM2**2.0)
    pan2.errorbar(time,tot,yerr=sigmaTot/np.sqrt(numPoints), fmt='ko', label=r'total')

    #pan1.plot(xrangefit*1e3,data_fitted,'g-')
    pan2.set_xlabel('rf freq')
    pan2.set_ylabel('counted number')

    fig3=plt.figure()
    pan3=fig3.add_subplot(1,1,1)
    pan3.errorbar(time,nPnorm,yerr=signPnorm/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
    pan3.errorbar(time,n0norm,yerr=sign0norm/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
    pan3.errorbar(time,nMnorm,yerr=signMnorm/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
    
    normTot= nMnorm+nPnorm+n0norm   
    sigmaNormTot=np.sqrt(signMnorm**2.0+sign0norm**2.0+signPnorm**2.0)
    if S==2:
        pan3.errorbar(time,nP2norm,yerr=signP2norm/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
        pan3.errorbar(time,nM2norm,yerr=signM2norm/np.sqrt(numPoints), fmt='mo', label=r'$m_F$=-2')
        normTot= nMnorm+nPnorm+n0norm+nP2norm+nM2norm   
        sigmaNormTot=np.sqrt(signMnorm**2.0+sign0norm**2.0+signPnorm**2.0+signP2norm*2.0+signM2norm**2.0)
    pan3.errorbar(time,normTot,yerr=sigmaNormTot/np.sqrt(numPoints), fmt='ko', label=r'total')

    #pan1.plot(xrangefit*1e3,data_fitted,'g-')
    pan3.set_xlabel('rf freq')
    pan3.set_ylabel('normalized counted number')
    
    fig4=plt.figure()
    pan4=fig4.add_subplot(1,1,1)
    pan4.errorbar(time,n1,yerr=sigma1/np.sqrt(numPoints), fmt='go',label='1st uwave shot')
    pan4.errorbar(time,n2,yerr=sigma2/np.sqrt(numPoints), fmt='bo',label='2nd uwave shot')
    pan4.errorbar(time,n1+n2,yerr=np.sqrt(sigma1**2.0+sigma2**2.0)/np.sqrt(numPoints), fmt='ko',label='total')
    pan4.legend()
    pan4.set_xlabel('rfFreq')
    pan4.set_ylabel('counted number')
   
    
    
    
    
fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 
pan1.errorbar(time,mag,yerr=sigmaMag,fmt='bo')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan1.set_xlabel('Time [s]')
pan1.set_ylabel('Magnetization')




fig2=plt.figure()
pan2=fig2.add_subplot(1,1,1)
pan2.errorbar(time*1e3,qlist,yerr=sigmaQ, fmt='bo')
pan2.set_xlabel('Oscillation time [ms]')
pan2.set_ylabel(r'Quasimomentum [$k_L$]')
stdGood = numPoints>1
(A,B,dA,dB)=lineFit.lineFit(time[:time.size]*1e3,qlist[:time.size],'Oscillation time [ms]',r'Quasimomentum [$k_L$]')#,errorbars=True,yerr=sigmaQ[stdGood],absolute_sigma=True,p0=(-0.26,0.78),maxfev=15000)
qlist=qlist-B
qlistFit=lineFit.line(time*1e3,A,0.0)

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(qlist,fractionP,yerr=sigmaP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(qlist,fraction0,yerr=sigma0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
pan1.errorbar(qlist,fractionM,yerr=sigmaM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan1.errorbar(qlist,fractionP2,yerr=sigmaP2/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
    pan1.errorbar(qlist,fractionM2,yerr=sigmaM2/np.sqrt(numPoints), fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan1.set_xlabel('quasimomentum [k_L]')
pan1.set_ylabel('Fractional populations')
legend()

fig1=plt.figure()
pan1=fig1.add_subplot(1,1,1)
#pan1.set_title(r'x0='+str(np.round(x0,6))+r', $\Gamma$='+str(np.round(Gamma,3))) 

pan1.errorbar(qlistFit,fractionP,yerr=sigmaP/np.sqrt(numPoints), fmt='bo', label=r'$m_F$=+1')
pan1.errorbar(qlistFit,fraction0,yerr=sigma0/np.sqrt(numPoints), fmt='go', label=r'$m_F$=0')
pan1.errorbar(qlistFit,fractionM,yerr=sigmaM/np.sqrt(numPoints), fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan1.errorbar(qlistFit,fractionP2,yerr=sigmaP2/np.sqrt(numPoints), fmt='co', label=r'$m_F$=+2')
    pan1.errorbar(qlistFit,fractionM2,yerr=sigmaM2/np.sqrt(numPoints), fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan1.set_xlabel('quasimomentum from fit [k_L]')
pan1.set_ylabel('Fractional populations')
legend()

if qlist.max()>0:
    newQlist=np.linspace(0.0,1.0,num=15)
else:
    newQlist=np.linspace(-1.0,0.0,num=15)
    
newP,newSigmaP=bestFitNewX(newQlist,qlistFit,fractionP,sigmaQ/np.sqrt(numPoints),sigmaP/np.sqrt(numPoints))
new0,newSigma0=bestFitNewX(newQlist,qlistFit,fraction0,sigmaQ/np.sqrt(numPoints),sigma0/np.sqrt(numPoints))
newM,newSigmaM=bestFitNewX(newQlist,qlistFit,fractionM,sigmaQ/np.sqrt(numPoints),sigmaM/np.sqrt(numPoints))
if S==2:
    newP2,newSigmaP2=bestFitNewX(newQlist,qlistFit,fractionP2,sigmaQ/np.sqrt(numPoints),sigmaP2/np.sqrt(numPoints))
    newM2,newSigmaM2=bestFitNewX(newQlist,qlistFit,fractionM2,sigmaQ/np.sqrt(numPoints),sigmaM2/np.sqrt(numPoints))

fig2=plt.figure()
pan=fig2.add_subplot(111)    
pan.errorbar(newQlist,newP,yerr=newSigmaP, fmt='bo', label=r'$m_F$=+1')
pan.errorbar(newQlist,new0,yerr=newSigma0, fmt='go', label=r'$m_F$=0')
pan.errorbar(newQlist,newM,yerr=newSigmaM, fmt='ro', label=r'$m_F$=-1')
if S==2:
    pan.errorbar(newQlist,newP2,yerr=newSigmaP2, fmt='co', label=r'$m_F$=+2')
    pan.errorbar(newQlist,newM2,yerr=newSigmaM2, fmt='mo', label=r'$m_F$=-2')
#pan1.plot(xrangefit*1e3,data_fitted,'g-')
pan.set_xlabel('new quasimomentum')
pan.set_ylabel('best fit fractional populations')
legend()    

np.savez(saveFileName,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,
         qlist=qlist,sigmaQ=sigmaQ,tlist=time,sigmaP=sigmaP,sigma0=sigma0,sigmaM=sigmaM,numPoints=numPoints,
         mag=mag,sigmaMag=sigmaMag,qlistFit=qlistFit,newQlist=newQlist,newP=newP,
         newSigmaP=newSigmaP,new0=new0,newSigma0=newSigma0,newM=newM,newSigmaM=newSigmaM)
if S==2:
    np.savez(saveFileName,fractionP2=fractionP2,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,fractionM2=fractionM2,
         qlist=qlist,sigmaQ=sigmaQ,tlist=time,sigmaP2=sigmaP2,sigmaP=sigmaP,sigma0=sigma0,sigmaM=sigmaM,sigmaM2=sigmaM2,numPoints=numPoints,
         mag=mag,sigmaMag=sigmaMag,qlistFit=qlistFit,newQlist=newQlist,newP=newP,
         newSigmaP=newSigmaP,new0=new0,newSigma0=newSigma0,newM=newM,newSigmaM=newSigmaM,
         newP2=newP2, newSigmaP2=newSigmaP2,newM2=newM2,newSigmaM2=newSigmaM2)
         
if getNums:
    np.savez(saveFileName,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,numP=nP,
             numM=nM, num0=n0, numPnorm=nPnorm, numMnorm=nMnorm,num0norm=n0norm,signM=signM,sign0=sign0,
             signP=signP,signMnorm=signMnorm,sign0norm=sign0norm,signPnorm=signPnorm,qlist=qlist,
             sigmaQ=sigmaQ,tlist=time,sigmaP=sigmaP,sigma0=sigma0,sigmaM=sigmaM,
             numPoints=numPoints,mag=mag,sigmaMag=sigmaMag,n1=n1,sigma1=sigma1,n2=n2,sigma2=sigma2,qlistFit=qlistFit)
    if S==2:
        np.savez(saveFileName,fractionP2=fractionP2,fractionP=fractionP,fraction0=fraction0,fractionM=fractionM,fractionM2=fractionM2,
                 numP2=nP2,numP=nP,numM=nM,numM2=nM2, num0=n0,numP2norm=nP2norm, numPnorm=nPnorm, 
                 numMnorm=nMnorm,numM2norm=nM2norm,num0norm=n0norm,signM2=signM2,signM=signM,sign0=sign0,signP=signP,signP2=signP2,
                 signM2norm=signM2norm,signMnorm=signMnorm,sign0norm=sign0norm,signPnorm=signPnorm,signP2norm=signP2norm,
                 qlist=qlist,sigmaQ=sigmaQ,tlist=time,sigmaP=sigmaP,sigmaP2=sigmaP2,sigma0=sigma0,sigmaM=sigmaM,numPoints=numPoints,
                 mag=mag,sigmaMag=sigmaMag,sigmaM2=sigmaM2,n1=n1,sigma1=sigma1,n2=n2,sigma2=sigma2,qlistFit=qlistFit)        
#fig3=plt.figure()
#pan3=fig3.add_subplot(1,1,1)
#pan3.errorbar(time*1e3, num1Tot,yerr=sigma1Tot,fmt= 'bo')
#pan3.errorbar(time*1e3, num2Tot,yerr=sigma2Tot,fmt= 'go')
#pan3.errorbar(time*1e3, num3Tot,yerr=sigma3Tot,fmt= 'ro')
#pan3.set_xlabel('Arp length [A]')
#pan3.set_ylabel('Total counted atom number')
#dataFile1=np.load('14Apr2016_files_254-298.npz')
#
#time1=dataFile1['tlist']
#sort1=np.argsort(time1)
#time1=time1[sort1]
#
#fractionP1=dataFile1['fractionP'][sort1]
#fraction01=dataFile1['fraction0'][sort1]
#fractionM1=dataFile1['fractionM'][sort1]
#qlist1=dataFile1['qlist'][sort1]
#for i in range(qlist1.size-1):
#    if qlist1[i+1]>qlist1[i]+0.5:
#        qlist1[i+1]=qlist1[i+1]-2.0
#        
#dataFile2=np.load('14Apr2016_files_299-343.npz')
#
#time2=dataFile2['tlist']
#sort2=np.argsort(time2)
#time2=time2[sort2]
#
#fractionP2=dataFile2['fractionP'][sort2]
#fraction02=dataFile2['fraction0'][sort2]
#fractionM2=dataFile2['fractionM'][sort2]
#qlist2=dataFile2['qlist'][sort2]
#for i in range(qlist2.size-1):
#    if qlist2[i+1]>qlist2[i]+0.5:
#        qlist2[i+1]=qlist2[i+1]-2.0
#        
#        
#dataFile3=np.load('14Apr2016_files_344-388.npz')
#
#time3=dataFile3['tlist']
#sort3=np.argsort(time3)
#time3=time3[sort3]
#
#fractionP3=dataFile3['fractionP'][sort3]
#fraction03=dataFile3['fraction0'][sort3]
#fractionM3=dataFile3['fractionM'][sort3]
#qlist3=dataFile3['qlist'][sort3]
#for i in range(qlist3.size-1):
#    if qlist3[i+1]>qlist3[i]+0.5:
#        qlist3[i+1]=qlist3[i+1]-2.0
#
#dataFile4=np.load('15Apr2016_files_129-173.npz')
#
#time4=dataFile4['tlist']
#sort4=np.argsort(time4)
#time4=time4[sort4]
#
#fractionP4=dataFile4['fractionP'][sort4]
#fraction04=dataFile4['fraction0'][sort4]
#fractionM4=dataFile4['fractionM'][sort4]
#qlist4=dataFile4['qlist'][sort4]
#for i in range(qlist4.size-1):
#    if qlist4[i+1]>qlist4[i]+0.5:
#        qlist4[i+1]=qlist4[i+1]-2.0
#        
#time=time1
#fractionP=(fractionP1+fractionP2+fractionP3+fractionP4)/4.0
#fraction0=(fraction01+fraction02+fraction03+fraction04)/4.0
#fractionM=(fractionM1+fractionM2+fractionM3+fractionM4)/4.0
#qlist=(qlist1+qlist2+qlist3+qlist4)/4.0

