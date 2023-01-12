#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 17:41:57 2021

@author: felipe
"""

import sys, os
sys.path.append(os.path.abspath('../spectrum'))
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pandas as pd
import Data
import scipy.signal as signal
from matplotlib import rc
rc('text', usetex=True)

def handle_close(evt):
        sys.exit()
class VData(object):
    def __init__(self,data,dataN,dataR,dataS,markers,norm,lims,fs):
        
        self.meanData=np.mean(data,axis=1)
        self.meanDataN=np.mean(dataN,axis=1)
        self.data=np.reshape(data,(np.shape(data)[0],16,16))
        self.dataN=np.reshape(dataN,(np.shape(dataN)[0],16,16))
        self.norm=norm
        self.minimo=lims[0]
        self.maximo=lims[1]
        self.fs=fs
        self.dataR=dataR
        self.dataS=dataS
        self.markers=markers
        self.fig = plt.figure(figsize=(8,6),frameon=False)
        gs1=gridspec.GridSpec(4,2,figure=self.fig,width_ratios=[0.3,0.7],height_ratios=[1,1,1,1],wspace=0.5)
        self.axim1=self.fig.add_subplot(gs1[0:2,0],frame_on=False)
        self.axim2=self.fig.add_subplot(gs1[2:4,0],frame_on=False)
        self.ax1=self.fig.add_subplot(gs1[0,1],frame_on=True)
        self.ax2=self.fig.add_subplot(gs1[1,1],frame_on=True)
        self.ax3=self.fig.add_subplot(gs1[2,1],frame_on=True)
        self.ax4=self.fig.add_subplot(gs1[3,1],frame_on=True)
        clb=plt.colorbar(cm.ScalarMappable(norm=self.norm, cmap=plt.cm.RdYlBu),ax=self.axim1,fraction=0.05)
        clb1=plt.colorbar(cm.ScalarMappable(norm=self.norm, cmap=plt.cm.RdYlBu),ax=self.axim2,fraction=0.05)
    def make_frame(self,t):
        n=int(t*self.fs)
        t_start=np.floor(t/6)
        n_second=int(t_start*6*self.fs)
        n_now=n-n_second
        tarray=np.linspace(0,6,6*self.fs+1)
        self.axim1.clear()
        self.axim2.clear()
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()
        self.axim1.tick_params(axis='both',which='both',length=0)
        self.axim1.imshow(self.data[n,:,:].T,cmap=plt.cm.RdYlBu,interpolation='bilinear',norm=self.norm)
        self.axim1.set_xticklabels([])
        self.axim1.set_yticklabels([])
        self.axim1.set_ylabel('$\phi_e$ z-score',fontsize=12)
        
        self.axim2.tick_params(axis='both',which='both',length=0)
        self.axim2.imshow(self.dataN[n,:,:].T,cmap=plt.cm.RdYlBu,interpolation='bilinear',norm=self.norm)
        self.axim2.set_xticklabels([])
        self.axim2.set_yticklabels([])
        self.axim2.set_ylabel('$\phi_n$ z-score',fontsize=12)
        
        self.ax1.plot(tarray[0:n_now],self.meanData[n_second:n_second+n_now]) 
        self.ax1.set_ylabel('$\phi_e$ z-score',fontsize=12)
        self.ax1.set_ylim([self.minimo,self.maximo])
        self.ax1.set_xlim([0,6])
        self.ax1.set_xticks(np.arange(0,8,2))
        self.ax1.set_xticklabels(np.arange(t_start*6,t_start*6+8,2))
        
        
        self.ax2.plot(tarray[0:n_now],self.dataR[n_second:n_second+n_now]) 
        self.ax2.set_ylabel('$\phi_r$ z-score',fontsize=12)
        self.ax2.set_ylim([self.minimo,self.maximo])
        self.ax2.set_xlim([0,6])
        self.ax2.set_xticks(np.arange(0,8,2))
        self.ax2.set_xticklabels(np.arange(t_start*6,t_start*6+8,2))
        
        self.ax3.plot(tarray[0:n_now],self.dataS[n_second:n_second+n_now]) 
        self.ax3.set_ylabel('$\phi_s$ z-score',fontsize=12)
        self.ax3.set_ylim([self.minimo,self.maximo])
        self.ax3.set_xlim([0,6])
        self.ax3.set_xticks(np.arange(0,8,2))
        self.ax3.set_xticklabels(np.arange(t_start*6,t_start*6+8,2))
        
        self.ax4.plot(tarray[0:n_now],self.meanDataN[n_second:n_second+n_now]) 
        self.ax4.plot(tarray[0:n_now],self.markers[n_second:n_second+n_now]+self.minimo) 
        self.ax4.set_ylabel('$\phi_n$ z-score',fontsize=12)
        self.ax4.set_ylim([self.minimo,self.maximo])
        self.ax4.set_xlim([0,6])
        self.ax4.set_xticks(np.arange(0,8,2))
        self.ax4.set_xticklabels(np.arange(t_start*6,t_start*6+8,2))
        self.ax4.set_xlabel('time (s)',fontsize=12)
        if n==0:
            plt.show()
        return mplfig_to_npimage(self.fig)
        
    
def main():
    if len(sys.argv) > 0:
        fs=100
        fps=60
        #Load txt file
        filename = sys.argv[1]
        filenameStim=re.sub('Phi','Stim',filename)
        filenameVideo=re.sub('Phi.txt','Video.mp4',filename)
        print('Loading and plotting '+filename+' ...')
        data, datazScore,std=Data.loadStoredData(filename)
        data, dataNzScore,std=Data.loadStoredData(filename,L=np.arange(256,512))
        data, dataSzScore,std=Data.loadStoredData(filename,L=[512])
        data, dataRzScore,std=Data.loadStoredData(filename,L=[513])
        dataStim, markers, time=Data.loadStoredStim(filenameStim,L=256,searchMarkers=True)
        data=datazScore.to_numpy(dtype=float)
        dataN=dataNzScore.to_numpy(dtype=float)
        dataR=dataRzScore.to_numpy(dtype=float)
        dataS=dataSzScore.to_numpy(dtype=float)
        lowFrequency=0.1
        highFrequency=40
        b,a=signal.cheby1(4,1e-6,[lowFrequency/(fs/2), highFrequency/(fs/2)],'bandpass')
        data=signal.filtfilt(b,a,data,axis=0)
        minimo=np.min(np.mean(data,axis=1))
        maximo=np.max(np.mean(data,axis=1))
        if minimo<-maximo:
            minimo=minimo//1.5*2
            maximo=-minimo
        else:
            maximo=maximo//1.5*2
            minimo=-maximo
        lims=(minimo,maximo)
        norm= colors.Normalize(vmin=minimo, vmax=maximo)
        a=VData(data,dataN,dataR,dataS,markers,norm,lims,fs)        
        time=np.shape(data)[0]
        print(time)
        anim = VideoClip(a.make_frame, duration=int((time//fs)))
        anim.write_videofile(filenameVideo, fps=fps)
        a.fig.canvas.mpl_connect('close_event',handle_close)
    else:
        print("Error: filename not valid")
if __name__=='__main__':
    main()        

