#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 18:38:11 2021

@author: felipe
"""

#imports
import sys, os
sys.path.append(os.path.abspath('../spectrum'))
sys.path.append(os.path.abspath('../plotting'))
sys.path.append(os.path.abspath('../eventsDetection'))
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import pandas as pd
import Data
import DataStrengths
import re
import Epochs
import spindles
import regression
import SO
import coincidences

#Process AAS simulations
def main():
    
    
    #Time , frequency Constants
    fs=100
    samples_coincidence=25
    h=1e-4
    h1=1e-2
    #Initialize arrays
    Spindles_coincidences_time=np.zeros((2000,))
    SOarray=np.zeros((2000,400))
    RMSarray=np.zeros((2000,400))
    SOextracted=0
    nsamples=0
    nstart=0
    flag_hold=False
    countC=0
    SpindlesSOcount=0
    flag_baseline=False
    if len(sys.argv) > 1:
        #Load data
        filename=sys.argv[1]
        print('Collecting Data from ',filename)
        filename_strengths=re.sub('Phi.txt','Strengths.txt',filename)
        tokensFilename=re.split('-',filename)
        if filename.find('baseline')!=-1:
            flag_baseline=True
            try:
        	    seed=int(tokensFilename[5])-1
            except:
                seed=0
        else:
            try:
                seed=int(tokensFilename[6])-1
            except:
                seed=0
        print(seed+1)
        Ephi=Data.loadStoredDataLong(filename,L=np.arange(0,256),maxLength=535000)
        if flag_baseline:
            meanBaseline=np.mean(Ephi)
            stdBaseline=np.std(Ephi)
            Ephiz=(Ephi-meanBaseline)/stdBaseline
        else:
            ##STORED BASELINE Thresholds events and statistics
            fileTh=np.load('BaselineEventsTh.npz')
            SOthreshold=fileTh['SOthreshold']
            SPthreshold=fileTh['SPthreshold']
            meanBaseline=fileTh['meanBaseline']
            stdBaseline=fileTh['stdBaseline']
            Ephiz=(Ephi-meanBaseline[seed])/stdBaseline[seed]
            
        data_length=len(Ephi)
        #Calculate loop-gains
        strengths, phi, ratios=DataStrengths.loadStoredStrengths(filename_strengths)
        G=DataStrengths.calculateGains(strengths,phi)
        X,Y,Z=DataStrengths.calculateXYZ(G)
        timeXYZ=np.arange(len(X))*2/60
        #Events with manual set baseline threshold 
        ##Change this, read from previous file
        
        # Baselines
        if flag_baseline:
            filteredSp, rmsSp,detectionSpindles, segmentationSpindles, countSpindles, SPthreshold=spindles.detectAPESpindlesBaseline(Ephiz,fs=fs,gap_min=0.1)
            detectionSO, SOsegmentation, SOs,countSO,SOthreshold=SO.detectAPESOBaseline(Ephi,neg_threshold=-1e-6,p2p_threshold=0)
            filteredD, rmsD,detectionDelta, countDelta=spindles.detectBand(Ephiz,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=SPthreshold[0])
        else:
            filteredSp, rmsSp,detectionSpindles, segmentationSpindles, countSpindles=spindles.detectAPESpindles(Ephiz,SPthreshold[:,seed],fs=fs,gap_min=0.1)
            detectionSO, SOsegmentation, SOs,countSO=SO.detectAPESO(Ephi,SOthreshold[:,seed],neg_threshold=-1e-6,p2p_threshold=0)
            filteredD, rmsD,detectionDelta, countDelta=spindles.detectBand(Ephiz,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=SPthreshold[0,seed])
        #hilbertDelta=signal.hilbert(filteredD)
        #phase=np.angle(hilbertDelta*np.exp(-1j*np.pi/2))+np.pi
        #Events per minute
        SPxminute,SPpercentage=spindles.eventsxMinute(segmentationSpindles)
        SOxminute,SOpercentage=spindles.eventsxMinute(SOsegmentation)
        #Coinicidences
        
        P_SO,P_SP,P_C=coincidences.eventsProbability(detectionSO,detectionSpindles,data_length=data_length,samplesCoincidence=25,sizeArray=2000)
        SpindlesSOcount, countC, Spindles_delay_time, coincidences_time=coincidences.coincidencesMeasures(detectionSpindles,detectionSO,SOsegmentation,segmentationSpindles,SOs,data_length=data_length,samplesCoincidence=25,h1=1e-2,sizeArray=2000)    
        for center in SOs:
            if center>200 and center<data_length-200:
                try:
                    SOarray[SOextracted,:]=filteredD[center-200:center+200]
                    RMSarray[SOextracted,:]=rmsSp[center-200:center+200,0]
                    SOextracted+=1
                except ValueError:
                    print(center)
        
        #Save
        del Ephi, filteredD,filteredSp,detectionSO
        freqs,scales,Sxx=Data.waveletComplexLong(Ephiz,maxLength=535000)
        print('Storing')
        if flag_baseline:
            np.savez(filename+'.npz',timeXYZ=timeXYZ,X=X,Y=Y,Z=Z,
                 SPxminute=SPxminute,SOxminute=SOxminute,
                 countSO=countSO,countSP=countSpindles,countC=SpindlesSOcount, countDelta=countDelta,
                 P_SO=P_SO,P_SP=P_SP,P_C=P_C,G=G,
                 Spindles_delay_time=Spindles_delay_time,
                 coincidences_time=coincidences_time,
                 strengths=strengths,
                 phi=phi,
                 ratios=ratios,
                 SOarray=SOarray,
                 SParray=RMSarray,
                 SOextracted=SOextracted,
                 freqs=freqs,
                 Sxx=Sxx,
                 SOthreshold=SOthreshold,
                 SPthreshold=SPthreshold,
                 rmsSp=rmsSp,
                 segmentationSpindles=segmentationSpindles,
                 meanBaseline=meanBaseline,
                 stdBaseline=stdBaseline
                 )
        else:
            np.savez(filename+'.npz',timeXYZ=timeXYZ,X=X,Y=Y,Z=Z,
                 SPxminute=SPxminute,SOxminute=SOxminute,
                 countSO=countSO,countSP=countSpindles,countC=SpindlesSOcount, countDelta=countDelta,
                 P_SO=P_SO,P_SP=P_SP,P_C=P_C,G=G,
                 Spindles_delay_time=Spindles_delay_time,
                 coincidences_time=coincidences_time,
                 strengths=strengths,
                 phi=phi,
                 ratios=ratios,
                 SOarray=SOarray,
                 SParray=RMSarray,
                 SOextracted=SOextracted,
                 freqs=freqs,
                 Sxx=Sxx,
                 rmsSp=rmsSp,
                 segmentationSpindles=segmentationSpindles
                 )
        del Ephiz, segmentationSpindles
        
if __name__=='__main__':
    main()
