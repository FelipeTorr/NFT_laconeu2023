#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu July 02 00:20:34 2021

@author: felipe
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import Wavelets
import sys

def loadStoredData1Ch(filename,L=np.arange(0,256),fs=100,sep='\t'):
	
	textfile = open(filename,'r')
	data=pd.read_csv(textfile, sep=sep,header=None,on_bad_lines='skip')
	data=data.to_numpy()
	data=np.array(data[:,L],dtype='float')
	meand=np.nanmean(data)
	cols=[]
	if len(L)>1:
		for j in L:
	    		cols.append('Node '+str(j+1))
		data1=np.copy(data)
		data=pd.DataFrame(data,columns=cols)
		data2=np.nanmean(data1,axis=1)
		stdd=np.nanstd(data2)
		try:
			datazScore=(data2-meand)/stdd
		except TypeError:
			datazScore=(data2-meand)
			print('Data type not allow z-scoring or standard deviation is zero')
		
	else:
		stdd=np.nanstd(data)
		if stdd>0:
			datazScore=(data-meand)/stdd
		else:
			datazScore=(data-meand)
	textfile.close()
	return data, datazScore,stdd

def loadStoredDataLong(filename,L=np.arange(0,256),fs=100,sep='\t',maxLength=541000):
    textfile = open(filename,'r')
    chunksize=5*fs
    nblocks=maxLength//(chunksize)
    data=np.array([])
    lowF=0.25
    highF=40
    b,a=signal.cheby1(4,1e-6,[lowF/(fs/2), highF/(fs/2)],'bandpass')
    df1 = pd.read_csv(textfile,sep=sep, header = None, skiprows=0, chunksize=chunksize,on_bad_lines='skip')
    for nindex in range(nblocks):
        df=df1.get_chunk()
        df=df.to_numpy()
        df=np.array(df[:,L],dtype='float')
        meand=np.nanmean(df,axis=1)
        data=np.hstack((data,meand))
    data=signal.filtfilt(b,a,data)
    textfile.close()
    return data
def loadStoredDataLongNoFilter(filename,L=np.arange(0,256),fs=100,sep='\t',maxLength=541000):
    textfile = open(filename,'r')
    chunksize=5*fs
    nblocks=maxLength//(chunksize)
    data=np.array([])
    df1 = pd.read_csv(textfile,sep=sep, header = None, skiprows=0, chunksize=chunksize,on_bad_lines='skip')
    for nindex in range(nblocks):
        df=df1.get_chunk()
        df=df.to_numpy()
        df=np.array(df[:,L],dtype='float')
        meand=np.nanmean(df,axis=1)
        data=np.hstack((data,meand))
    textfile.close()
    return data

def loadStoredData(filename,L=np.arange(0,256),fs=100):
	
	textfile = open(filename,'r')
	data=pd.read_csv(textfile, sep='\t')
	data=data.to_numpy()
    
	data=np.array(data[0::,L],dtype='float')
	cols=[]
	for j in L:
		cols.append('Node '+str(j+1))
	data=pd.DataFrame(data)
	std=stats.tstd(data.to_numpy(dtype=float))
	data1=data.apply(signal.detrend)
	try:
		datazScore=data1.apply(stats.zscore)
	except TypeError:
		datazScore=data1;
		print('Data type not allow z-scoring or standard deviation is zero')
	textfile.close()
	
	return data, datazScore,std
	
def loadStoredStim(filename,L=256,searchMarkers=False,storageCorrection=False):
	textfile = open(filename,'r')
	data=pd.read_csv(textfile, sep='\t')
	data=data.to_numpy()
	if storageCorrection:
		lengthdata=np.shape(data)[0]//200*200-1
		for i in np.arange(lengthdata,399,-200):
			data[i,:]=data[i-200,:]
	marker=np.array(data[:,L])#,dtype='int8')
	meandata=np.mean(data[:,0:L-2],axis=1)
	if searchMarkers:
		marker=np.zeros_like(marker)
		hold=0
		for n in range(np.shape(data)[0]):
			if meandata[n] != 0 and hold==0:
                #perturbation detected
				marker[n]=1
				#wait until perturbation end
				hold=1
			elif meandata[n] == 0 and hold==1:
				#ready for next perturbation
				hold=0;
	
	time=np.array(data[:,L+1],dtype='float')
	#Warning: data from now retain only nodes data
	data=data[:,0:L]
	textfile.close()

	return meandata, marker, time

def loadStoredStimLong(filename,L=256,searchMarkers=False,storageCorrection=False,maxLength=541000):
    textfile = open(filename,'r')
    chunksize=500
    nblocks=maxLength//(chunksize)
    data=np.array([])
    marker=np.array([])
    time=np.array([])
    sep='\t'
    df1 = pd.read_csv(textfile,sep=sep, header = None, skiprows=0, chunksize=chunksize,on_bad_lines='skip')
    for nindex in range(nblocks):
        df=df1.get_chunk()
        df=df.to_numpy()
        ddata=np.array(df[:,0:L],dtype='float')
        dmarker=np.array(df[:,L])
        dtime=np.array(df[:,L+1],dtype='float')
        meand=np.nanmean(ddata,axis=1)
        data=np.hstack((data,meand))
        marker=np.hstack((marker,dmarker))
        time=np.hstack((time,dtime))
        
    if storageCorrection:
        lengthdata=np.shape(data)[0]//200*200-1
        for i in np.arange(lengthdata,399,-200):
            data[i,:]=data[i-200,:]
    meandata=data
    if searchMarkers:
        marker=np.zeros_like(marker)
        hold=0
        for n in range(np.shape(data)[0]):
            if meandata[n] != 0 and hold==0:
                #perturbation detected
                marker[n]=1
				#wait until perturbation end
                hold=1
            elif meandata[n] == 0 and hold==1:
				#ready for next perturbation
                hold=0;
    textfile.close()

    return meandata, marker, time
    

def loadStoredPhase(filename,fs=100,sep=','):
	
    textfile = open(filename,'r')
    data=pd.read_csv(textfile, sep=sep,header=None,on_bad_lines='skip')
    data=data.to_numpy()
    data=np.array(data,dtype='float')
    onlineFiltered=data[:,0] 
    onlineEnvelope=data[:,1] 
    onlinePhase=data[:,2]
        
    	
    return onlineFiltered, onlineEnvelope, onlinePhase

def loadStoredPhaseLong(filename,fs=100,sep=',',maxLength=541000):
    filtered=np.array([])
    envelope=np.array([])
    phase=np.array([])
    chunksize=500
    nblocks=maxLength//(chunksize)
    textfile = open(filename,'r')
    df1 = pd.read_csv(textfile,sep=sep, header = None, skiprows=0, chunksize=chunksize,on_bad_lines='skip')
    for nindex in range(nblocks):
        df=df1.get_chunk()
        df=df.to_numpy()
        dfiltered=np.array(df[:,0],dtype='float')
        denvelope=np.array(df[:,1],dtype='float')
        dphase=np.array(df[:,2],dtype='float')
        filtered=np.hstack((filtered,dfiltered))
        envelope=np.hstack((envelope,denvelope))
        phase=np.hstack((phase,dphase))

    return filtered, envelope, phase

def spatialFilter(data,sigma=0.141421,Nx=16,Ny=16):
	lengthTime=np.shape(data)[0]
	filtered=np.zeros((Nx,Ny,lengthTime))
	for n in range(lengthTime):
		data_n=np.transpose(np.reshape(data[n,:],(Nx,Ny))) #Transpose to keep the same order of Matlab reshape
		filtered[:,:,n]=ndimage.gaussian_filter(data_n,sigma=sigma,mode='wrap')
	return filtered
	
def meanSpatialFilter(data,sigma=0.141421,Nx=16,Ny=16):
	lengthTime=np.shape(data)[0]
	filtered=spatialFilter(data,sigma=sigma,Nx=Nx,Ny=Ny)
	return np.transpose(np.mean(np.reshape(filtered,(Nx*Ny,lengthTime)),axis=0))

def loadBaseline(filenameBase,order=True,filter='Uniform',fs=100):
	#Load SHAM data
	textBase = open(filenameBase,'r')
	dataB=pd.read_csv(textBase, sep='\t')
	dataB=dataB.to_numpy()
	LBE=np.arange(0,256)
	LBN=np.arange(256,512)
	dataBE=np.array(dataB[:,LBE],dtype='float')
	dataBN=np.array(dataB[:,LBN],dtype='float')
	
    #Spatial mean or filtering
	if filter=='Uniform':
		dataBE=np.mean(dataBE,axis=1)
		dataBN=np.mean(dataBN,axis=1)
	elif filter=='Gaussian':
		dataBE=meanSpatialFilter(dataBE,sigma=0.141)
		dataBN=meanSpatialFilter(dataBN,sigma=0.141)
        
    #Baseline statistics
	meanE=np.mean(dataBE)
	meanN=np.mean(dataBN)
	stdE=stats.tstd(dataBE)
	stdN=stats.tstd(dataBN)
	dataBE=dataBE-meanE
	dataBN=dataBN-meanN
	dataBEz=dataBE/stdE
	dataBNz=dataBN/stdN
    
	return dataBE,dataBN,dataBEz,dataBNz

def loadAllData(filename,filenameBase,order=True,filter='Uniform',fs=100):
	
	filenameStim=filename[0:-7]+'Stim.txt'
	
	#Load STIM data
	textfile = open(filename,'r')
	data=pd.read_csv(textfile, sep='\t')
	data=data.to_numpy()
	textfile.close()
    #Load SHAM data
	textBase = open(filenameBase,'r')
	dataB=pd.read_csv(textBase, sep='\t')
	dataB=dataB.to_numpy()
    #Load stimulation
	dataS, marker, time=loadStoredStim(filenameStim,searchMarkers=True)
	textBase.close()
    
    #Separate Phi_e from Phi_n
	LE=np.arange(0,256)
	LN=np.arange(256,512)
	LBE=np.arange(0,256)
	LBN=np.arange(256,512)
	if order==True:
		dataE=np.array(data[:,LE],dtype='float')
		dataN=np.array(data[:,LN],dtype='float')
	else:
		dataE=np.array(data[:,LBE],dtype='float')
		dataN=np.array(data[:,LBN],dtype='float')
	
    #Shapes-Again	
	#dataBE=np.array(dataB[:,LE],dtype='float')
	#dataBN=np.array(dataB[:,LN],dtype='float')
	dataBE=np.array(dataB[:,LBE],dtype='float')
	dataBN=np.array(dataB[:,LBN],dtype='float')
	
	#Spatial mean or filtering
	if filter=='Uniform':
		dataE=np.nanmean(dataE,axis=1)
		dataN=np.nanmean(dataN,axis=1)
		dataBE=np.nanmean(dataBE,axis=1)
		dataBN=np.nanmean(dataBN,axis=1)
	elif filter=='Gaussian':
		dataE=meanSpatialFilter(dataE,sigma=0.141)
		dataN=meanSpatialFilter(dataN,sigma=0.141)
		dataBE=meanSpatialFilter(dataBE,sigma=0.141)
		dataBN=meanSpatialFilter(dataBN,sigma=0.141)
	

	#Baseline statistics
	meanE=np.nanmean(dataBE)
	meanN=np.nanmean(dataBN)
	stdE=np.nanstd(dataBE)
	stdN=np.nanstd(dataBN)
	#Detrend
	dataE=dataE-meanE
	dataN=dataN-meanN
	dataBE=dataBE-meanE
	dataBN=dataBN-meanN
	
	#z-Score (divsion by global std)
	dataEz=dataE/stdE
	dataN=dataN/stdN
	dataBE=dataBE/stdE
	dataBN=dataBN/stdN
	
	return dataE, dataEz, dataN, dataBE, dataBN, dataS, marker, time



def spectrogram(datazScore,fs=100,correctF=False):
	window=signal.get_window('hamming', int(10*fs))
	noverlap=int(9.5*fs)
	meanSignal=np.nanmean(datazScore.to_numpy(dtype=float),axis=1)
	f,t,Sxx=signal.spectrogram(meanSignal, fs=fs, window=window, noverlap=noverlap)
	if correctF:
		f=f[1::]
		for col in range(np.shape(Sxx)[1]):
			Sxx[1::,col]=Sxx[1::,col]*f
	return f,t,Sxx

def spectrogramSingle(meanSignal,fs=100,correctF=False):
	window=signal.get_window('hamming', int(10*fs))
	noverlap=int(9.5*fs)
	f,t,Sxx=signal.spectrogram(meanSignal, fs=fs, window=window, noverlap=noverlap)
	if correctF:
		f=f[1::]
		for col in range(np.shape(Sxx)[1]):
			Sxx[1::,col]=Sxx[1::,col]*f
	return f,t,Sxx

def wavelet(datazScore,fs=100, correctF=False):
	meanSignal=np.nanmean(datazScore.to_numpy(dtype=float),axis=1)
# 	mother = pycwt.Morlet(6)
# 	s0 = 2/fs  # Starting scale
# 	octaves_suboctaves = 1 / 8 # Twelve sub-octaves per octaves
# 	N = 7 / octaves_suboctaves  # Eigth powers of two per total sub-octave
# 	coefs, scales, freqs, coi, fft, fftfreqs = pycwt.cwt(meanSignal,1/fs,octaves_suboctaves,s0,N,mother)
	freqs = np.logspace(-1,1.478,300)
# 	freqs=np.arange(0.1,30.1,0.1)
	dt=1/fs
	wavel=Wavelets.Morlet(meanSignal,freqs=freqs,dt=dt,omega0=15)
    
	coefs=wavel.getnormpower()
	scales=wavel.getscales()
	if correctF:
		for col in range(np.shape(coefs)[1]):
			coefs[:,col]=coefs[:,col]/scales[:]
	return np.flip(freqs),np.flip(scales),np.flipud(coefs)
# 	return freqs,scales,coefs

def waveletSingle(meanSignal,fs=100, correctF=False):
	freqs = np.logspace(-1,1.478,300)
	dt=1/fs
	wavel=Wavelets.Morlet(meanSignal,freqs=freqs,dt=dt,omega0=15)
    
	coefs=wavel.getnormpower()
	scales=wavel.getscales()
	if correctF:
		for col in range(np.shape(coefs)[1]):
			coefs[:,col]=coefs[:,col]*freqs[:]
	return np.flip(freqs),np.flip(scales),np.flipud(coefs)


def waveletComplex(datazScore,fs=100, correctF=False):
	meanSignal=np.nanmean(datazScore.to_numpy(dtype=float),axis=1)
	#mother = pycwt.Morlet(6)
	#s0 = 2/fs  # Starting scale
	#octaves_suboctaves = 1 / 8 # Twelve sub-octaves per octaves
	#N = 7 / octaves_suboctaves  # Eigth powers of two per total sub-octave
	#coefs, scales, freqs, coi, fft, fftfreqs = pycwt.cwt(meanSignal,1/fs,octaves_suboctaves,s0,N,mother)
	freqs = np.logspace(-1,1.478,300)

	dt=1/fs
	wavel=Wavelets.Morlet(meanSignal,freqs=freqs,dt=dt,omega0=15)
    
	coefs=wavel.getdata()
	scales=wavel.getscales()	

	if correctF:
		for col in range(np.shape(coefs)[1]):
			coefs[:,col]=coefs[:,col]/freqs[:]
	return np.flip(freqs),np.flip(scales),np.flipud(coefs)
# 	return freqs,scales,coefs

def waveletComplexSingle(meanSignal,fs=100, correctF=False):
	freqs = np.logspace(-1,1.478,300)
	dt=1/fs
    
	wavel=Wavelets.Morlet(meanSignal,freqs=freqs,dt=dt,omega0=15)
    
	coefs=wavel.getdata()
	scales=wavel.getscales()	

	if correctF:
		for col in range(np.shape(coefs)[1]):
			coefs[:,col]=coefs[:,col]/freqs[:]
	return np.flip(freqs),np.flip(scales),np.flipud(coefs)

def waveletComplexLong(meanSignal,fs=100, maxLength=541000,correctF=False):
    freqs = np.logspace(-1,1.478,300)
    dt=1/fs
    chunkSize=10*fs
    
    nblocks=maxLength//chunkSize    
    spectrums=np.zeros((nblocks,300))
	
    for nchunk in range(nblocks-1):
        wavel=Wavelets.Morlet(meanSignal[nchunk*chunkSize:(nchunk+1)*chunkSize],freqs=freqs,dt=dt,omega0=15)
        coefs=wavel.getdata()
        if nchunk==1:
            scales=wavel.getscales()
        if correctF:
            for col in range(np.shape(coefs)[1]):
                coefs[:,col]=coefs[:,col]/freqs[:]
        spectrums[nchunk,:]=np.flipud(np.sum(np.abs(coefs),axis=1))
    return np.flip(freqs),np.flip(scales),spectrums

def welch(datazScore,fs=100,startTime=0,endTime=0):
	window=signal.get_window('hamming', int(10*fs))
	noverlap=int(8*fs)
	meanSignal=np.nanmean(datazScore.to_numpy(dtype=float),axis=1)
	if endTime==0:
		f,Pxx=signal.welch(meanSignal, fs=fs, window=window, noverlap=noverlap)
	elif endTime>startTime:
		f,Pxx=signal.welch(meanSignal[(startTime*fs)//1:(endTime*fs)//1+1], fs=fs, window=window, noverlap=noverlap)
	else:
		print('EndTime must be higher than starTime, lower than the signal length o zero to get all time periodogram')
		f=0
		Pxx=-1
	return f, Pxx

def welchSingle(meanSignal,fs=100,startTime=0,endTime=0):
	window=signal.get_window('hamming', int(10*fs))
	noverlap=int(8*fs)
	if endTime==0:
		f,Pxx=signal.welch(meanSignal, fs=fs, window=window, noverlap=noverlap)
	elif endTime>startTime:
		f,Pxx=signal.welch(meanSignal[(startTime*fs)//1:(endTime*fs)//1+1], fs=fs, window=window, noverlap=noverlap)
	else:
		print('EndTime must be higher than starTime, lower than the signal length o zero to get all time periodogram')
		f=0
		Pxx=-1
	return f, Pxx
	

def main():
	if len(sys.argv) > 0:
		fs=100
		filename = sys.argv[1]
		data, datazScore=loadStoredData(filename)
		meanSignal=np.nanmean(data.to_numpy(dtype=float),axis=1)
		meanSignalz=np.nanmean(datazScore.to_numpy(dtype=float),axis=1)
		freqs,scales,Sxx,coi=wavelet(datazScore,correctF=False)
		freqs1,t,Sxx1=spectrogram(datazScore)
		filtered=spatialFilter(data)
		filteredZ=spatialFilter(datazScore)
		print(freqs)
		print(np.shape(t)[0])
		plt.figure()
		plt.plot(meanSignal)
		plt.plot(filtered[0,0,:])
		plt.figure()
		plt.plot(meanSignalz)
		plt.plot(np.mean(filteredZ,axis=(0,1)))
		plt.figure()
		plt.imshow(np.abs(Sxx),aspect='auto')
		plt.figure()
		plt.imshow(10*np.log10(Sxx1),aspect='auto')
		plt.show()
		
		
	else:
		print ("Error: you must input the filename")
		sys.exit()
if __name__=='__main__':
	main()
	

	
	
