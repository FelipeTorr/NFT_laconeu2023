#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:50:11 2021

@author: felipe
"""
import numpy as np

meansSOp2p=np.array([[0.0010868276440743384,0.0010878097188067972],
           [0.004085888499127412,0.039357973413892126],
           [0.008272821436812473,0.008272821436812473]])

rmsSPs=np.array([[0.7891359580443994,0.8374597155898004],
        [0.17911885360887597,0.22698121051622014],
        [0.30091519696511404,0.30091519696511404]])


np.savez('BaselinesPlasticity.npz',meansSOp2p=meansSOp2p,rmsSPs=rmsSPs)