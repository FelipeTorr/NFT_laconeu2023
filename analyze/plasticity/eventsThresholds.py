import numpy as np

SOthreshold=np.zeros((541000,5))
SPthreshold=np.zeros((541000,5))
meanBaseline=np.zeros((5,))
stdBaseline=np.zeros((5,))
for seed in range(1,6):
    file=np.load('/home/felipe/Dropbox/UTFSM/NeuralField/analyze/collectedData/APE/SWSAPE1-Phin1-1.5H-2pi3-seed-%d-baseline-Phi.txt.npz'%seed)
    SOthreshold[:,seed-1:seed]=file['SOthreshold']
    SPthreshold[:,seed-1:seed]=file['SPthreshold']
    meanBaseline[seed-1]=file['meanBaseline']
    stdBaseline[seed-1]=file['stdBaseline']
#SPth=np.array([0.6111836521778009,0.6099439810047651,0.6143083995836686,0.6110427786797491,0.6174907686225877])
#SOth=np.array([0.01606903090572213,0.01589077216838762,0.01576979837878471,0.01601856311349591,0.015927606282061104])


np.savez('BaselineEventsTh.npz', SPthreshold=SPthreshold,SOthreshold=SOthreshold,meanBaseline=meanBaseline,stdBaseline=stdBaseline)

