#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 11:16:47 2021

@author: felipe
"""

import numpy as np

def detectCoincidences(detectionSpindles,detectionSO,samplesCoincidence=25,sizeArray=2000):
    coincidences_time=np.zeros((sizeArray,))
    tcoincidences=detectionSO+detectionSpindles[:,0]
    coincidence=np.zeros((len(tcoincidences)))
    nsamples=0
    nstart=0
    flag_hold=False
    jj=0
    samples_coincidence=samplesCoincidence
    for ii in range(len(tcoincidences)):
        if tcoincidences[ii]==2:
            nsamples+=1
            if flag_hold==False:
                nstart=ii
                flag_hold=True
        if tcoincidences[ii]<2:
            if flag_hold:
                flag_hold=False
                if nsamples>=samples_coincidence:
                    # print(nstart,ii,(ii-nstart)*0.01)
                    coincidence[nstart:ii]=1
                    nsamples=0
                    coincidences_time[jj]=nstart
                    jj+=1
            else:
                nsamples=0
    coincidence=np.argwhere(coincidence==1)
    return coincidence, coincidences_time
    
def coincidencesMeasures(detectionSpindles,detectionSO,segmentationSO,segmentationSP,SOs,data_length=9000,samplesCoincidence=25,h1=1e-2,sizeArray=2000):
    dSp=detectionSpindles[1:len(detectionSpindles),0]-detectionSpindles[0:-1,0]
    upB=np.argwhere(dSp==1)[:,0]
    downB=np.argwhere(dSp==-1)[:,0]
    centersB=(upB+downB)//2
    Spindles_delay_time=np.zeros((sizeArray,))
    SpindlesSOcount=0
    countC=0
    coincidence, coincidences_time=detectCoincidences(detectionSpindles, detectionSO,samplesCoincidence=samplesCoincidence,sizeArray=sizeArray)
    
    if len(coincidence)>0:   
        m=0
        mm=0
        u=-1 
        for coinc in coincidence:
            if segmentationSP[coinc]!=u:
                u=int(segmentationSP[coinc][0][0])-1
                Spindles_delay_time[SpindlesSOcount+mm:SpindlesSOcount+mm+1]=(SOs[int(segmentationSO[coinc][0])-1]-centersB[u])*h1
                mm+=1
                u+=1
            if coinc in list(centersB):
                m+=1
        SpindlesSOcount+=mm
        countC+=m
    return SpindlesSOcount, countC, Spindles_delay_time, coincidences_time*h1

def eventsProbability(detectionSO,detectionSpindles,data_length=9000,samplesCoincidence=25,sizeArray=2000):

    coincidence, coincidences_time=detectCoincidences(detectionSpindles, detectionSO,samplesCoincidence=samplesCoincidence,sizeArray=sizeArray)
    P_SO=len(np.nonzero(detectionSO==1)[0])/data_length
    P_SP=len(np.nonzero(detectionSpindles==1)[0])/data_length
    P_C=len(coincidence)/data_length
    
    return P_SO,P_SP,P_C
        