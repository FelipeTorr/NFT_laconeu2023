#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:58:30 2021

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
    ##Thresholds events
    file=np.load('BaselinesPlasticity.npz')
    meansSOp2p=file['meansSOp2p']
    rmsSPs=file['rmsSPs']
    #Time , frequency Constants
    fs=100
    samples_coincidence=25
    h=1e-4
    h1=1e-2
    lowF=0.1
    highF=40
    b,a=signal.cheby1(4,1e-6,[lowF/(fs/2), highF/(fs/2)],'bandpass')
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
    if len(sys.argv) > 1:
        #Load data
        filename=sys.argv[1]
        print('Collecting Data from ',filename)
        filename_strengths=re.sub('Phi.txt','Strengths.txt',filename)
        Ephi, Ez,stdE=Data.loadStoredData1Ch(filename,L=np.arange(0,256))
        data_length=np.shape(Ez)[0]
        Ephi=Ephi.to_numpy()
        Ephi=np.mean(Ephi,axis=1)-np.mean(Ephi)
        #Detrend data
        #detrendPolyEphi,trendPolyEphi=regression.detrend(Ephi,samples=500,iterations=10)
        detrendPolyEphi=signal.filtfilt(b,a,Ephi)
        detrendPolyEz=(detrendPolyEphi-np.mean(detrendPolyEphi))/np.std(detrendPolyEphi)
        #Calculate loop-gains
        strengths, phi, ratios=DataStrengths.loadStoredStrengths(filename_strengths)
        G=DataStrengths.calculateGains(strengths,phi)
        X,Y,Z=DataStrengths.calculateXYZ(G)
        timeXYZ=np.arange(len(X))*2/60
        
        plt.plot(detrendPolyEphi)
        
        #Events with manual set baseline threshold 
        ##Change this, read from previous file
        ## Baselines
        if filename.find('baseline')!=-1:
            detectionSO, SOsegmentation, SOs,countSO,meanSOp2p=SO.SODetection(detrendPolyEphi,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=0)
            filteredSp, rmsSp,detectionSpindles, segmentationSpindles, countSpindles=spindles.detectsegmentSpindles(detrendPolyEz,fs=fs,threshold=0)
            meanRMSSp=np.mean(rmsSp)
            print(meanSOp2p)
            print(np.mean(rmsSp))
        else:
            phin=int(filename[filename.find('Phin')+4])-1
            #seed=int(filename[filename.find('seed')+5])-1
            if filename.find('N2')!=-1:
                
                meanSOp2p=meansSOp2p[0,phin]
                meanRMSSp=rmsSPs[0,phin]
            elif filename.find('N3')!=-1:
                meanSOp2p=meansSOp2p[1,phin]
                meanRMSSp=rmsSPs[1,phin]
            elif filename.find('SWS')!=-1:
                meanSOp2p=meansSOp2p[2,phin]
                meanRMSSp=rmsSPs[2,phin]
        ###
        meanRMSD=1e-2
        meanRMST=1e-2
        
        detectionSO, SOsegmentation, SOs,countSO,meanSOp2p1=SO.SODetection(detrendPolyEphi,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=1.25*meanSOp2p)
        filteredSp, rmsSp,detectionSpindles, segmentationSpindles, countSpindles=spindles.detectsegmentSpindles(detrendPolyEz,fs=fs,threshold=1.25*meanRMSSp)
        filteredD, rmsD,detectionDelta, countDelta=spindles.detectBand(detrendPolyEz,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=1.25*meanRMSD)
        filteredT, rmsT,detectionTheta, countTheta =spindles.detectBand(detrendPolyEz,fs=fs,lowFrequency=4.5,highFrequency=8,threshold=1.25*meanRMST)
        
        hilbertDelta=signal.hilbert(filteredD)
        phase=np.angle(hilbertDelta*np.exp(-1j*np.pi/2))+np.pi
        #Events per minute
        SPxminute,SPpercentage=spindles.eventsxMinute(segmentationSpindles)
        SOxminute,SOpercentage=spindles.eventsxMinute(SOsegmentation)
        #Coinicidences
        plt.plot(0.015*detectionSO+0.06)
        plt.plot(0.015*detectionSpindles+0.08)
        plt.savefig('TimeDomain.pdf',dpi=300)
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
        print('Storing')
        #Save
        np.savez(filename+'.npz',timeXYZ=timeXYZ,X=X,Y=Y,Z=Z,
                 SPxminute=SPxminute,SOxminute=SOxminute,
                 countSO=countSO,countSP=countSpindles,countC=SpindlesSOcount, countDelta=countDelta,
                 countTheta=countTheta,
                 P_SO=P_SO,P_SP=P_SP,P_C=P_C,G=G,
                 Spindles_delay_time=Spindles_delay_time,
                 coincidences_time=coincidences_time,
                 strengths=strengths,
                 phi=phi,
                 ratios=ratios,
                 SOarray=SOarray,
                 SParray=RMSarray,
                 SOextracted=SOextracted,
                 meanRMSSp=meanRMSSp,
                 meanRMSD=rmsD,
                 meanRMST=rmsT,
                 meanSOp2p=meanSOp2p)
if __name__=='__main__':
    main()
