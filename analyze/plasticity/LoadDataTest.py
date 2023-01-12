#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 02:19:50 2021

@author: felipe
"""

import sys, os
sys.path.append(os.path.abspath('../spectrum'))
sys.path.append(os.path.abspath('../eventsDetection'))
sys.path.append(os.path.abspath('../plotting'))
import Data
import SO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as signal


# filename1='/media/felipe/Elements/APESWS/SWSAPE1-Phin1-1.5H-2pi3-seed-3-baseline-Phi.txt'
filename2='/media/felipe/Elements/APESWS/openloop/SWSAPE1-Phin1-StimTimes2-1.5H-2pi3-seed-3-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-Poisson-Phi.txt'
# filename2='/media/felipe/Elements/APESWS/openloop/SWSAPE1-Phin1-StimTimes2-1.5H-2pi3-seed-1-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-Poisson-Phi.txt'
# filename3='/media/felipe/Elements/APESWS/openloop/SWSAPE1-Phin1-StimTimes3-1.5H-2pi3-seed-3-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-OFF-Phi.txt'
#'/media/felipe/Elements/APESWS/openloop/SWSAPE1-Phin1-StimTimes3-1.5H-2pi3-seed-1-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-OFF-Phi.txt'
# filename4='/media/felipe/Elements/APESWS/openloop/SWSAPE1-Phin1-StimTimes3-1.5H-2pi3-seed-3-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-OFF-Phi.txt'
# phie1=Data.loadStoredDataLong(filename1,L=np.arange(0,256),maxLength=540000)
phie2=Data.loadStoredDataLong(filename2,L=np.arange(0,256),maxLength=540000)
# phie3=Data.loadStoredDataLong(filename3,L=np.arange(0,256),maxLength=5000)
# phie4=Data.loadStoredDataLong(filename4,L=np.arange(0,256),maxLength=20000)
#%%
# plt.plot(phie1[300000:400000])
plt.figure(figsize=(7.5,3.5))
plt.plot(phie2[0:2000])
plt.plot(segmentationSpindles[0:2000]*0.001)
plt.ylim([-0.015,0.03])
# plt.plot(phie3[1000:2000])

# plt.plot(phie4[0:200])
# phiez=(phie1-np.mean(phie1))/np.std(phie1)
    # plt.plot(phie)
    # freqs,scales,Sxx=Data.waveletComplexSingle(phiez)
# freqs,scales,SxxL=Data.waveletComplexLong(phiez,maxLength=7500)
# plt.figure()
    # plt.imshow(np.abs(SxxL[:,33:216].T),aspect='auto')
    # plt.figure()
    # plt.loglog(freqs[33:216],np.sum(np.abs(Sxx[33:216,:]),axis=1))
# plt.loglog(freqs[33:216],np.sum(np.abs(SxxL[:,33:216].T),axis=1))
    # h1=1e-2
    
    
    # print('Events detection')
    # detectionSO,SOs,countSO=SO.SODetection(phie)
    # b,a=signal.cheby1(4,0.001,[0.5/(fs/2), 1.25/(fs/2)],'bandpass')
    # filtered_D=signal.filtfilt(b,a,Phiez)
    # stdD=np.std(filtered_D)
    # del filtered_D
    # b,a=signal.cheby1(4,0.001,[4.5/(fs/2), 8.0/(fs/2)],'bandpass')
    # filtered_T=signal.filtfilt(b,a,Phiez)
    # stdT=np.std(filtered_T)
    # del filtered_T
    # b,a=signal.cheby1(4,0.001,[11/(fs/2), 16/(fs/2)],'bandpass')
    # filtered_S=signal.filtfilt(b,a,Phiez)
    # stdS=np.std(filtered_S)
    # del filtered_S

    # detectionSO,SOs,countSO=SO.SODetection(Phie)
    # filtered, rms,detectionSpindles, countSpindles=spindles.detectSpindles(Phiez,fs=fs)  
    # filtered, rms,detectionTheta, countTheta=spindles.detectBand(Phiez,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=.26,tmin=0.45,gapMax=0.05,windowLength=0.1)
    # filtered, rms,detectionDelta, countDelta=spindles.detectBand(Phiez,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=.26,tmin=1.0,gapMax=0.05,windowLength=0.1)
    
