#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 00:44:42 2021

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

#Process AAS simulations
def main():
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
    if len(sys.argv) > 1:
        #Load data
        filename=sys.argv[1]
        print('Collecting Data from ',filename)
        filename_stim=re.sub('Phi.txt','Stim.txt',filename)
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
        Ephi=Data.loadStoredDataLong(filename,L=np.arange(0,256),maxLength=480000)
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
        data, stimMarkers, time=Data.loadStoredStim(filename_stim,L=256,searchMarkers=True,storageCorrection=True)
        #Phase same online filter
        del data, time
        offlineFiltered=signal.filtfilt(b1,a1,Ephiz,axis=0)
        offlinePhase=np.angle(signal.hilbert(offlineFiltered)*np.exp(-1j*np.pi/2))+np.pi
        markers=np.nonzero(stimMarkers)
        phaseStim=offlinePhase[markers]
        markers=markers[0]
        markersNumber=len(markers)
        print(markersNumber)
        for stim in range(markersNumber):
            pos_marker=markers[stim]
            if phaseStim[stim]>-7/4*np.pi and phaseStim[stim]<=np.pi/4:
                phases0.append((pos_marker,phaseStim[stim]))
                data0=np.hstack((data0,Ephiz[pos_marker-200:pos_marker+600]))
            elif phaseStim[stim]>1/4*np.pi and phaseStim[stim]<=np.pi*3/4:
                phases90.append((pos_marker,phaseStim[stim]))
                data90=np.hstack((data90,Ephiz[pos_marker-200:pos_marker+600]))
            elif phaseStim[stim]>3/4*np.pi and phaseStim[stim]<=np.pi*5/4:
                phases180.append((pos_marker,phaseStim[stim]))
                data180=np.hstack((data180,Ephiz[pos_marker-200:pos_marker+600]))
            elif phaseStim[stim]>5/4*np.pi and phaseStim[stim]<=np.pi*7/4:
                phases270.append((pos_marker,phaseStim[stim]))
                data270=np.hstack((data270,Ephiz[pos_marker-200:pos_marker+600]))
        phases0=np.array(phases0)
        phases90=np.array(phases90)
        phases180=np.array(phases180)
        phases270=np.array(phases270)
        print('finished:',filename)
        print(len(phases0),len(phases90),len(phases180),len(phases270))
        np.savez(filename+'-phaseStim.npz',phases0=phases0,phases90=phases90,
                 phases180=phases180,phases270=phases270,
                 data0=data0,
                 data90=data90,data180=data180,data270=data270)

if __name__=='__main__':
    main()