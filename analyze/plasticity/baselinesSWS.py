#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 16:12:17 2021

@author: felipe
"""

#imports
import sys, os
sys.path.append(os.path.abspath('../spectrum'))
sys.path.append(os.path.abspath('../plotting'))
import DataStrengths
import numpy as np

import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import rc
import steadyFunctions
rc('text', usetex=True)

fig=plt.figure(figsize=(7.5,7.5))
gs1=gridspec.GridSpec(4, 2,hspace=0.3,wspace=0.4)
axA=fig.add_subplot(gs1[0,0])
axB=fig.add_subplot(gs1[1,0])
axC=fig.add_subplot(gs1[2,0])
axD=fig.add_subplot(gs1[3,0])
axE=fig.add_subplot(gs1[0,1])
axF=fig.add_subplot(gs1[1,1])
axG=fig.add_subplot(gs1[2,1])
axH=fig.add_subplot(gs1[3,1])

fig2=plt.figure(figsize=(7.5,7.5))
gs2=gridspec.GridSpec(4, 2,hspace=0.3,wspace=0.4)
ax2A=fig2.add_subplot(gs2[0,0])
ax2B=fig2.add_subplot(gs2[1,0])
ax2C=fig2.add_subplot(gs2[2,0])
ax2D=fig2.add_subplot(gs2[3,0])
ax2E=fig2.add_subplot(gs2[0,1])
ax2F=fig2.add_subplot(gs2[1,1])
ax2G=fig2.add_subplot(gs2[2,1])
ax2H=fig2.add_subplot(gs2[3,1])

fig3=plt.figure(figsize=(7.5,7.5))
gs3=gridspec.GridSpec(5, 1,hspace=0.4,wspace=0.4)
ax3A=fig3.add_subplot(gs3[0])
ax3B=fig3.add_subplot(gs3[1])
ax3C=fig3.add_subplot(gs3[2])
ax3D=fig3.add_subplot(gs3[3])
ax3E=fig3.add_subplot(gs3[4])

fig4=plt.figure(figsize=(7.5,7.5))
gs4=gridspec.GridSpec(4, 2,hspace=0.4,wspace=0.5)
ax4A=fig4.add_subplot(gs4[0,0])
ax4B=fig4.add_subplot(gs4[1,0])
ax4C=fig4.add_subplot(gs4[2,0])
ax4D=fig4.add_subplot(gs4[3,0])
ax4E=fig4.add_subplot(gs4[0,1])
ax4F=fig4.add_subplot(gs4[1,1])
ax4G=fig4.add_subplot(gs4[2,1])
ax4H=fig4.add_subplot(gs4[3,1])

fig5=plt.figure(figsize=(8,7.5))
gs5=gridspec.GridSpec(1, 3,wspace=0.4,hspace=0.4)
axA5=fig5.add_subplot(gs5[0])
axB5=fig5.add_subplot(gs5[1])
axC5=fig5.add_subplot(gs5[2])

fig6=plt.figure(figsize=(7.5,5.2))
gs6=gridspec.GridSpec(2,2,hspace=0.4,wspace=0.5)
axA6=fig6.add_subplot(gs6[0,0])
axB6=fig6.add_subplot(gs6[0,1])
axC6=fig6.add_subplot(gs6[1,0])
axD6=fig6.add_subplot(gs6[1,1])

fig7=plt.figure(figsize=(7.5,7.5))
gs7=gridspec.GridSpec(2,2,hspace=0.4,wspace=0.5)
axA7=fig7.add_subplot(gs7[0,0],projection='3d')



#interpolation phin=0
phi0=np.zeros((101,5))

#interpolation phin=0.5
phi05=np.zeros((101,5))

#interpolation phin=1
phi1=np.zeros((101,5))

points=np.arange(101)

#Model parameters
Qmax=340
sigma_p=0.0038
nus=np.array([5.541, -5.652, 1.53, 0.286, 1.12, 2.67, -1.73, 9.22])*1e-3
nu1=np.array([3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3])
nu2=np.array([6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.6e-3])

alpha=45
beta=186
gamma=116
r=0.086

print("Load fsolve data")
# points_low=np.arange(21)
# for phinn in np.arange(0.25,2.7,0.25):
#     for step in range(1,102):
       
#         #h5py requires the absolute path
#         matfile=h5py.File('/home/felipe/Dropbox/UTFSM/NeuralField/analyze/collectedData/NF-interpolation/interpolation_step_'+str(step)+'-phin1.mat','r')
#         solutions_varios=np.array(matfile['soulutions_e'])
#         nusi=(np.array(matfile['nus'])).T
#         if phinn==0.25:
#             #Phin=0
#             index_phin=np.argwhere(solutions_varios[:,0]>0)[0][0]
#             phi0[step-1,:]=[solutions_varios[index_phin,1],solutions_varios[index_phin,1],solutions_varios[index_phin,2],solutions_varios[index_phin,3],solutions_varios[index_phin,0]]
            
#             #Phin=0
#             index_phin1=np.argwhere(solutions_varios[:,0]>=1)[0][0]
#             phi1[step-1,:]=[solutions_varios[index_phin1,1],solutions_varios[index_phin1,1],solutions_varios[index_phin1,2],solutions_varios[index_phin1,3],solutions_varios[index_phin1,0]]
    
#         #Phin=0.5
#         index_phin05=np.argwhere(solutions_varios[:,0]>=phinn)[0][0]
#         phi05[step-1,:]=[solutions_varios[index_phin05,1],solutions_varios[index_phin05,1],solutions_varios[index_phin05,2],solutions_varios[index_phin05,3],solutions_varios[index_phin05,0]]
#     G05=DataStrengths.calculateGains(nusi,phi05)
#     X05,Y05,Z05=DataStrengths.calculateXYZ(G05)
#     print(phinn,Y05[::5])
#     axA7.scatter(X05[::5], Y05[::5], Z05[::5], s=10,marker='o',c=points_low,cmap=cm.summer_r)
# # Poweri[:,step-1]=steadyFunctions.SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,nusi[:,step-1],t0,k0) 

# G0=DataStrengths.calculateGains(nusi,phi0)
# X0,Y0,Z0=DataStrengths.calculateXYZ(G0)

# G1=DataStrengths.calculateGains(nusi,phi1)
# X1,Y1,Z1=DataStrengths.calculateXYZ(G1)

# XX=np.arange(0.98,1.0,0.005)
# YY=1-XX
# ZZ=0
# axA7.scatter(XX, YY, np.zeros_like(XX), s=10,marker='s',c='k')

filenames=['/media/felipe/TOSHIBA2T/1h30m/SWSAAS-Phin1-1.5H-2pi3-seed-1-baseline-Strengths.txt',
           '/media/felipe/TOSHIBA2T/1h30m/SWSAAS-Phin1-1.5H-2pi3-seed-2-baseline-Strengths.txt',
           '/media/felipe/TOSHIBA2T/1h30m/SWSAAS-Phin1-1.5H-2pi3-seed-3-baseline-Strengths.txt',
           '/media/felipe/TOSHIBA2T/1h30m/SWSAAS-Phin1-1.5H-2pi3-seed-4-baseline-Strengths.txt',
           '/media/felipe/TOSHIBA2T/1h30m/SWSAAS-Phin1-1.5H-2pi3-seed-5-baseline-Strengths.txt']


casesStage=['SHAM1','SHAM2','SHAM3','SHAM4','SHAM5']
legendsStimulation=['SHAM','STIM-R','STIM-P']
for m,file in enumerate(filenames):
        strengths, phi, ratios=DataStrengths.loadStoredStrengths(file)
        G=DataStrengths.calculateGains(strengths,phi)
        X,Y,Z=DataStrengths.calculateXYZ(G)
        timeXYZ=np.arange(len(X))*2/60
        axA.plot(timeXYZ[0:200],ratios[0:200,0])
        axB.plot(timeXYZ[0:200],ratios[0:200,1])
        axC.plot(timeXYZ[0:200],ratios[0:200,2])
        axD.plot(timeXYZ[0:200],ratios[0:200,3])
        axE.plot(timeXYZ[0:200],ratios[0:200,4])
        axF.plot(timeXYZ[0:200],ratios[0:200,5])
        axG.plot(timeXYZ[0:200],ratios[0:200,6])
        axH.plot(timeXYZ[0:200],ratios[0:200,7])
            
        ax2A.plot(timeXYZ[0:200],strengths[0:200,0])
        ax2B.plot(timeXYZ[0:200],strengths[0:200,1])
        ax2C.plot(timeXYZ[0:200],strengths[0:200,2])
        ax2D.plot(timeXYZ[0:200],strengths[0:200,3])
        ax2E.plot(timeXYZ[0:200],strengths[0:200,4])
        ax2F.plot(timeXYZ[0:200],strengths[0:200,5])
        ax2G.plot(timeXYZ[0:200],strengths[0:200,6])
        ax2H.plot(timeXYZ[0:200],strengths[0:200,7])

        
        ax3A.plot(timeXYZ[0:200],phi[0:200,0])
        ax3B.plot(timeXYZ[0:200],phi[0:200,1])
        ax3C.plot(timeXYZ[0:200],phi[0:200,2])
        ax3D.plot(timeXYZ[0:200],phi[0:200,3])
        ax3E.plot(timeXYZ[0:200],phi[0:200,4],':')
        
        ax4A.plot(timeXYZ[0:200],G[0:200,0])
        ax4B.plot(timeXYZ[0:200],G[0:200,1])
        ax4C.plot(timeXYZ[0:200],G[0:200,2])
        ax4D.plot(timeXYZ[0:200],G[0:200,3])
        ax4E.plot(timeXYZ[0:200],G[0:200,4])
        ax4F.plot(timeXYZ[0:200],G[0:200,5])
        ax4G.plot(timeXYZ[0:200],G[0:200,6])
        ax4H.plot(timeXYZ[0:200],G[0:200,7])
        
        axA5.plot(timeXYZ[0:200],G[0:200,0],color=plt.cm.tab10(m),label=r'$G_{ee}$')
        axA5.plot(timeXYZ[0:200],G[0:200,1],'--',color=plt.cm.tab10(m),label=r'$G_{ei}$')
        axB5.plot(timeXYZ[0:200],G[0:200,2],color=plt.cm.tab10(m),label=r'$G_{es}$')
        axB5.plot(timeXYZ[0:200],G[0:200,5],'--',color=plt.cm.tab10(m),label=r'$G_{se}$')
        axB5.plot(timeXYZ[0:200],G[0:200,6],':',color=plt.cm.tab10(m),label=r'$G_{sr}$')
        axC5.plot(timeXYZ[0:200],G[0:200,7],color=plt.cm.tab10(m),label=r'$G_{sn}$')
        axC5.plot(timeXYZ[0:200],G[0:200,3],'--',color=plt.cm.tab10(m),label=r'$G_{re}$')
        axC5.plot(timeXYZ[0:200],G[0:200,4],':',color=plt.cm.tab10(m),label=r'$G_{rs}$')
        
        
        axA6.plot(timeXYZ[0:200],X[0:200])
        #axA6.plot(X[0],Y[0],'o',markersize=0.1,color=plt.cm.tab10(3*m))
        
        axB6.plot(timeXYZ[0:200],Y[0:200])
        #axB6.plot(X[0],Z[0],'o',markersize=0.1,color=plt.cm.tab10(3*m))
        
        axC6.plot(timeXYZ[0:200],Z[0:200])
        
        axD6.plot(1-(X[0:200,50]+Y[0:200,50]),Z[0:200,50])
        
        if m==4:
            axA7.scatter(X[0],Y[0],Z[0],marker='P',s=70,color='k')
            axA7.scatter(X[900],Y[900],Z[900],marker='P',s=60,color='b')
            axA7.scatter(X[1800],Y[1800],Z[1800],marker='P',s=60,color='r')
            axA7.scatter(X[0:],Y[0:],Z[0:],marker='o',s=10)
        

        
        

axA.set_ylabel(r'$s_{ee}$')
axB.set_ylabel(r'$s_{ei}$')
axC.set_ylabel(r'$s_{es}$')
axD.set_ylabel(r'$s_{re}$')
axE.set_ylabel(r'$s_{rs}$')
axF.set_ylabel(r'$s_{se}$')
axG.set_ylabel(r'$s_{sr}$')
axH.set_ylabel(r'$s_{sn}$')
axA.legend(casesStage,fontsize=8)    
ax2A.set_ylabel(r'$\nu_{ee}$')
ax2B.set_ylabel(r'$\nu_{ei}$')
ax2C.set_ylabel(r'$\nu_{es}$')
ax2D.set_ylabel(r'$\nu_{re}$')
ax2E.set_ylabel(r'$\nu_{rs}$')
ax2F.set_ylabel(r'$\nu_{se}$')
ax2G.set_ylabel(r'$\nu_{sr}$')
ax2H.set_ylabel(r'$\nu_{sn}$')
ax2H.legend(casesStage,ncol=2,fontsize=6,columnspacing=0.3)


ax3A.set_ylabel(r'$\phi_{e}$')
ax3B.set_ylabel(r'$\phi_{i}$')
ax3C.set_ylabel(r'$\phi_{r}$')
ax3D.set_ylabel(r'$\phi_{s}$')
ax3E.set_ylabel(r'$\phi_{n}$')
ax3E.set_xlabel('time (min.)')
ax3A.legend(casesStage,ncol=2,fontsize=6,columnspacing=0.3)

ax4A.set_ylabel('Gee')
ax4B.set_ylabel('Gei')
ax4C.set_ylabel('Ges')
ax4D.set_ylabel('Gre')
ax4E.set_ylabel('Grs')
ax4F.set_ylabel('Gse')
ax4G.set_ylabel('Gsr')
ax4H.set_ylabel('Gsn')
ax4H.legend(casesStage,ncol=2,fontsize=6,columnspacing=0.3)

axA5.set_xlabel('Time (min.)')
axB5.set_xlabel('Time (min.)')
axC5.set_xlabel('Time (min.)')

axA5.set_ylabel('$G_{ab}(t)-G_{ab}(t_0)$')
axB5.set_ylabel('$G_{ab}(t)-G_{ab}(t_0)$')
axC5.set_ylabel('$G_{ab}(t)-G_{ab}(t_0)$')
axA5.legend()
axB5.legend()
axC5.legend()

axA6.set_xlabel('Time (min.)')
axA6.set_ylabel('X')

axB6.set_xlabel('Time (min.)')
axB6.set_ylabel('Y')

axC6.set_xlabel('Time (min.)')
axC6.set_ylabel('Z')

axD6.set_xlabel('1-(X+Y)')
axD6.set_ylabel('Z')
axD6.legend(casesStage,fontsize=6,ncol=2,columnspacing=0.3)


#Interpolation
# axA7.scatter(X0, Y0, Z0, s=10,c=points,cmap=cm.summer_r,marker='o')
# axA7.scatter(X1, Y1, Z1, s=10,marker='o',c=points,cmap=cm.summer_r)
# axA7.text(X0[0], Y0[0]+0.002,Z0[0], r'N2 $\phi_n^{(0)}=0$', None)
# axA7.text(X0[100], Y0[100]+0.002, Z0[100], 'N3', None)
# # axA7.text(X05[0], Y05[0]+0.002,Z05[0], r'N2 $\phi_n^{(0)}=0.5$', None)
# # axA7.text(X05[100], Y05[100]+0.002, Z05[100], 'N3', None)
# axA7.text(X1[0], Y1[0]+0.002,Z1[0], r'N2 $\phi_n^{(0)}=1$', None)
# axA7.text(X1[100], Y1[100]+0.002, Z1[100], 'N3', None)

# axA7.scatter(X1[66], Y1[66], Z1[66], marker='d',s=70, c='black',edgecolors='black')


axA7.set_xlabel('X')
axA7.set_ylabel('Y')
axA7.set_zlabel('Z')

# axA7.view_init(30,0)
# fig7.savefig('BaselinesALL0.pdf',dpi=300,bbox_inches='tight')
# axA7.view_init(30,45)
# fig7.savefig('BaselinesALL1.pdf',dpi=300,bbox_inches='tight')
# axA7.view_init(30,90)
# fig7.savefig('BaselinesALL2.pdf',dpi=300,bbox_inches='tight')
# axA7.view_init(90,0)
# fig7.savefig('BaselinesALL3.pdf',dpi=300,bbox_inches='tight')


print(G[0,:])
print(X[0],Y[0],Z[0])

# fig3.savefig('BaselinesALLPhi.pdf',dpi=30,bbox_inches='tight')
# fig4.savefig('BaselinesPresolverGains.pdf',dpi=30,bbox_inches='tight')
# fig2.savefig('BaselinesPresolverNus.pdf',dpi=300,bbox_inches='tight')
#fig7.savefig('N2Phin1Baseline.pdf',dpi=300,bbox_inches='tight')

plt.show()