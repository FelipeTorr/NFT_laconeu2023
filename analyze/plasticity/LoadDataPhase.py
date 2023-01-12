#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 01:33:35 2021

@author: felipe
"""

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


filename='/home/felipe/timeseries/SWSAPE1-Phin1-StimTimes1-1.5H-2pi3-seed-1-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-ON-Phi.txt'
Ephi=Data.loadStoredDataLong(filename,L=np.arange(0,256),maxLength=30000)
lowF=0.82
highF=0.88
fs=100
b,a=signal.cheby1(4,1e-6,[lowF/(fs/2), highF/(fs/2)],'bandpass')
EphiSO=signal.filtfilt(b,a,Ephi)

filename='/home/felipe/timeseries/SWSAPE1-Phin1-StimTimes1-1.5H-2pi3-seed-1-p-decreaseRamp-F0.85-A3.46e+01-D5.00e-02-sdA-OFF-sdF-ON-Phase.txt'
Filtered=Data.loadStoredDataLongNoFilter(filename,L=[0],maxLength=30000)
Envelope=Data.loadStoredDataLongNoFilter(filename,L=[1],maxLength=30000)
Phase=Data.loadStoredDataLongNoFilter(filename,L=[2],maxLength=30000)
#%%
#plt.plot(Ephi[0:30000])
plt.plot(Filtered[8000:9000])
plt.plot(Envelope[8000:9000])
plt.plot(Phase[8000:9000]*0.005)