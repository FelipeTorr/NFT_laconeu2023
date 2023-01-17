#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 17:58:08 2021

@author: felipe
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

t=np.arange(-0.3,0.3,0.001)
etas=[10,50]
for n in range(2):
    H=np.zeros((len(t),2))
    A=-1
    Am=1
    eta=etas[n]
    
    H[np.argwhere(t<=0)[:,0],0]=Am*1.1*np.exp(eta*t[t<=0])
    H[np.argwhere(t<=0)[:,0],1]=A*np.exp(eta*t[t<=0])
    H[np.argwhere(t>0)[:,0],0]=Am*1.1*np.exp(-eta*t[t>0])
    H[np.argwhere(t>0)[:,0],1]=Am*np.exp(-eta*t[t>0])
    plt.subplot(1,2,n+1)
    plt.plot(t,H)
    plt.hlines(0,xmin=-0.3,xmax=0.3,color='k' )
    plt.legend(['CDP','STDP' ])
    plt.xlabel('tiempo (s)')
    plt.ylabel(r'h($\eta t$)')
    plt.title(r'$\eta$=%d'%eta)
plt.tight_layout()
plt.savefig('plasticidad.pdf',dpi=300)