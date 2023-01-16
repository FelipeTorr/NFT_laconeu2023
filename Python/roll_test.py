#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:57:33 2019

@author: felipe
"""

import numpy as np
import time
import matplotlib.pyplot as plt
N=100
length=1000
short_history=np.zeros((100,length))
long_history=np.zeros((100,length*5))
result1=np.zeros((100,))
result2=np.zeros((100,))
samples_delay=500
stored1=[]
stored2=[]

storing_index=[]
delay_index=[]
#Set m at the m-step, carry one position
now=time.time()
for j in range(5):
    for m in range(length):
        short_history[:,1]=m+j*1000
        result1+=0.1*short_history[:,samples_delay]
        stored1.append(result1[0])
        short_history=np.roll(short_history,1,axis=1)
        
elapsed1=time.time()-now

#Set m%length at the m-step position
now1=time.time()
real_index=0;
samples_delay_now=0
for m in range(5*length):
    new_value=m%length
    j=m//length
    store_index=(m+samples_delay)%(5*length) #Guardar tiene un slot menos
    long_history[:,store_index]=new_value+j*1000
    result2+=0.1*long_history[:,samples_delay_now];  
    stored2.append(result2[0])
    storing_index.append(store_index)
    delay_index.append(samples_delay_now)
    samples_delay_now+=1
    if samples_delay_now>=5*length:
        samples_delay_now=0
elapsed2=time.time()-now1    

plt.figure()
plt.plot(np.array(storing_index)-np.array(delay_index))
plt.show()