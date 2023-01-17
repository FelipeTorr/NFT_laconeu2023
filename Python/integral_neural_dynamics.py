#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 12:11:22 2023

@author: felipe
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

#%%
#Time-frequency parameters
fs=1000
duration=0.2
t=np.linspace(0,0.2,int(duration*fs//1))

#Sleep
alpha=45
beta=186
L1=alpha*beta/(beta-alpha)*(np.exp(-alpha*t)-np.exp(-beta*t))

#Awake
alpha=83
beta=769
L2=alpha*beta/(beta-alpha)*(np.exp(-alpha*t)-np.exp(-beta*t))

plt.plot(t,L1,'k',label='Sleep')
plt.plot(t,L2,label='Awake')
plt.xlabel('time (s)')
plt.ylabel('Dendritic response function (s^{-1})')
plt.legend()
plt.savefig('dendritic_response.pdf',dpi=300)

#%%
#Let convolve
phi0=np.zeros_like(t)/fs;
phi0[20]=1/fs;
phi0[100]=1/fs;
phi0[180]=1/fs;

phi1=np.ones_like(t)/fs;

V0=np.convolve(L1, phi0)
nu0=1e-3

V1=np.convolve(L1, phi1)
nu1=1e-3

t1=np.linspace(0,2*duration-1/fs,int((2*duration-1/fs)*fs))
plt.figure()
plt.subplot(3,1,1)
plt.plot(t,phi0*0.1*fs,label='0.1 input')
plt.plot(t1,V0,label='output')
plt.xlabel('time (s)')
plt.ylabel('D1 (s^{-1})')
plt.legend()
plt.subplot(3,1,2)
plt.plot(t,phi1*fs,label='input')
plt.plot(t1,V1,label='output')
plt.xlabel('time (s)')
plt.ylabel('D2 (s^{-1})')
plt.legend()
plt.subplot(3,1,3)
#Soma voltage
somaV=nu0*V0+nu1*V1
plt.plot(somaV,label='Soma voltage')
plt.xlabel('time (s)')
plt.ylabel('Soma Voltage (V)')
plt.legend()
plt.tight_layout()
plt.savefig('soma_voltage.pdf',dpi=300)

#%%
#Sigmoid
v=np.arange(-0.02,0.04,0.001)
def sigmoid(V,Qmax=340,theta=0.01292,sigma_p=0.0038):
    firing_response=Qmax/(1+np.exp(-(V-theta)/sigma_p))
    return firing_response

def inverseSigmoid(Q,Qmax=340,theta=0.01292,sigma_p=0.0038):
    iQ=theta-sigma_p*(np.log(Qmax-Q)-np.log(Q))
    return iQ
#errf
def errf(V,Qmax=340,theta=0.01292,sigma_p=0.0038):
    sigma=sigma_p*np.pi/np.sqrt(3)
    normal=stats.norm(loc=theta,scale=sigma)
    er=Qmax*normal.cdf(V)
    return er

def derivateSigmoid(V,Qmax=340,theta=0.01292,sigma_p=0.0038):
    derivate=sigmoid(V,Qmax=Qmax,theta=theta,sigma_p=sigma_p)/sigma_p*(1-sigmoid(V,Qmax=Qmax,theta=theta,sigma_p=sigma_p)/Qmax)
    return derivate

Q=sigmoid(v)
Q1=errf(v)

plt.figure()
plt.plot(v,Q,label='Sigmoid')
plt.plot(v,Q1,':.',label='errf')
plt.xlabel('Voltage (V)')
plt.ylabel('firing response (s^{-1})')


##linearization
v1=0.0
l1=sigmoid(v1)+derivateSigmoid(v1)*(v-v1)
plt.plot(v,l1,'--',label='Linear approx. at V=0.00')
##linearization
v2=0.01
l2=sigmoid(v2)+derivateSigmoid(v2)*(v-v2)
plt.plot(v,l2,'--',label='Linear approx. at V=0.01')

plt.legend()
plt.ylim([0,350])
plt.savefig('FiringResponse.pdf',dpi=300)
