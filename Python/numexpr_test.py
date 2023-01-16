#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:59:05 2019

@author: felipe
"""

import numpy as np
import numexpr as ne
import time

def multiplix(X,Y):
    M = X.shape[0]
    N = Y.shape[1]
    D = np.empty((M, N), dtype=np.float)
    for i in range(M):
        a=X[i,:]
        for j in range(N):
            b=Y[:,j]
            D[i,j]=ne.evaluate('sum(a*b)')
    return D
def sigmoid(x):
    Q=1
    s=5
    t=400
    y=ne.evaluate('Q/(1+exp(-(x-t)/s))')
    return y

def sigmoid1(x):
    Q=1
    s=5
    t=400
    y=Q/(1+np.exp(-(x-t)/s))
    return y

a=np.ones((800,4000))
b=np.ones((4000,2))
b[2,0]=5.45;
b[1,1]=8.12;
a[59,31]=12.24
a[54,20]=123.56899


#now=time.time()
#for m in range(10000):
#    c=multiplix(a,b)
#    y=sigmoid(c)
#elapsed=time.time()-now

now1=time.time()
for m in range(100):
    d=np.matmul(a,b)
    #w=sigmoid1(d)
elapsed1=time.time()-now1

#print(elapsed)
print(elapsed1)

##Conclusion numpy (in fact python standard kernel) is optimized for non-square matrix multiplication. Choose numpy for calculations.