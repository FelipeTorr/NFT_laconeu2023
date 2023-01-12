#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 21:08:23 2021

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
import statistic_plots as sp

#Process AAS baseline simulations

#Time , frequency Constants
fs=100
samples_coincidence=25
h=1e-4
h1=1e-2
#Filter
a1=np.array([1,-4.980302427753561,9.927029769710018,-9.899189875772134,4.938508551069478,-0.986045936040990,0])
b1=np.array([5.925095844680851e-10,1.135074946542897e-09,1.065592361719296e-10,-8.452568314723746e-10,-5.657541641346955e-10,-2.898181150705222e-10,-1.333146565053192e-10])
    
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
phases0=[]
phases90=[]
phases180=[]
phases270=[]
data0=np.array(([]))
data90=np.array(([]))
data180=np.array(([]))
data270=np.array(([]))
for seed in range(1,6):
    #Load data
    filename='/media/felipe/Elements/APESWS/baseline/SWSAPE1-Phin1-1.5H-2pi3-seed-%d-baseline-Phi.txt'%seed
    print('Collecting Data from ',filename)
    print(seed)
    Ephi=Data.loadStoredDataLong(filename,L=np.arange(0,256),maxLength=541000)
    meanBaseline=np.mean(Ephi)
    stdBaseline=np.std(Ephi)
    Ephiz=(Ephi-meanBaseline)/stdBaseline
    offlineFiltered=signal.filtfilt(b1,a1,Ephiz,axis=0)
    offlinePhase=np.angle(signal.hilbert(offlineFiltered)*np.exp(-1j*np.pi/2))+np.pi
    np.savez('Baseline-seed-%d-PhizPhase.npz'%seed,Ephiz=Ephiz,offlinePhase=offlinePhase)
    del offlineFiltered,offlinePhase,Ephi,Ephiz