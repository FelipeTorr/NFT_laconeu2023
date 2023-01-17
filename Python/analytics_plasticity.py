#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 22:38:42 2020

@author: felipe
"""

import numpy as np
from scipy.optimize import fsolve

def sigmoid(x,Qmax=340,theta=0.01292,sigma_p=0.0038):
    Q=Qmax/(1+np.exp((theta-x)/sigma_p))
    return Q



def derivateSigmoid(x,Qmax=340,theta=0.01292,sigma_p=0.0038):
    rho=sigmoid(x,Qmax=Qmax,theta=theta,sigma_p=sigma_p)*1/(sigma_p)*(1-sigmoid(x,Qmax=Qmax,theta=theta,sigma_p=sigma_p)*1/Qmax)
    return rho

def derivateSigmoidfiring(q,Qmax=340,theta=0.01292,sigma_p=0.0038):
    rho=q*1/(sigma_p)*(1-q*1/Qmax)
    return rho

def inverseSigmoid(x,Qmax=340,theta=0.01292,sigma_p=0.0038):
    try:
        voltage=np.zeros_like(x,dtype=float)
        for n in range(len(x)):
            if x[n]>0:
                voltage[n]=theta-sigma_p*(np.log(Qmax-x[n])-np.log(x[n]))
            else:
                voltage[n]=-np.exp(-x[n]) 
                #It must go to -infinity (exponentially is an option)
    except TypeError:
        if x>0:
            voltage=theta-sigma_p*(np.log(Qmax-x)-np.log(x))
        else:
            voltage=-np.exp(-x) 
    return voltage

def Dspatial(omega, gamma, r):
    k=0
    D=(1-1j*omega/gamma)**2+k**2*r**2
    return D

def Dtemporal(omega,alpha,beta):
    D=(1-1j*omega/alpha)*(1-1j*omega/beta);
    return D

    
def CharacteristicEquation(omega,k,alpha,beta,gamma,r,rho,strengths,t0):
    Gsrs=rho[2]*rho[1]*strengths[6]*strengths[4];
    Gee=rho[0]*strengths[0];
    Gei=rho[0]*strengths[1];
    Gese=rho[0]*rho[2]*strengths[2]*strengths[5];
    Gesre=rho[0]*rho[1]*rho[2]*strengths[2]*strengths[6]*strengths[3];
    L=1/Dtemporal(omega,alpha,beta)
    q2r2=(1-1j*omega/gamma)**2-(1/(1-Gei*L))*((Gee*L+((Gese*L**2+Gesre*L**3)*np.exp(1j*omega*t0))/(1-Gsrs*L**2)))
    s=(1-Gei*L)*(1-Gsrs*L**2)*(q2r2+k**2*r**2)
    return s

def steadystate(phie,strengths,phin=1):
    vee=strengths[0]
    vei=strengths[1]
    ves=strengths[2]
    vre=strengths[3]
    vrs=strengths[4]
    vse=strengths[5]
    vsr=strengths[6]
    vsn=strengths[7]
    #Other populations solutions
    phi_s=(1/ves)*inverseSigmoid(phie)-(vee+vei)*phie
    phi_r=sigmoid(vre*phie+vrs*phi_s)
    phi_sp=sigmoid(vse*phie+vsr*phi_r+vsn*phin)
    v=inverseSigmoid(phie)-(vee+vei)*phie-ves*phi_sp
    return v

def calculateRho(phin,strengths):
    phi_e=fsolve(steadystate,[10],(strengths,phin),xtol=1e-8)
    vee=strengths[0]
    vei=strengths[1]
    ves=strengths[2]
    vre=strengths[3]
    vrs=strengths[4]
    vse=strengths[5]
    vsr=strengths[6]
    vsn=strengths[7] 
    phi_sp=(1/ves)*inverseSigmoid(phi_e)-(vee+vei)*phi_e
    phi_r=sigmoid(vre*phi_e+vrs*phi_sp)
    phi_s=sigmoid(vse*phi_e+vsr*phi_r+vsn*phin)
    rho_e=derivateSigmoidfiring(phi_e,Qmax=340,theta=0.01292,sigma_p=0.0038)
    rho_s=derivateSigmoidfiring(phi_s,Qmax=340,theta=0.01292,sigma_p=0.0038)
    rho_r=derivateSigmoidfiring(phi_r,Qmax=340,theta=0.01292,sigma_p=0.0038)
    print(phi_e,phi_r,phi_s)
    rho=np.array([rho_e,rho_r,rho_s])
    return rho
    
def TransferFunction(omega,k,alpha,beta,gamma,r,rho,strengths,t0):
    Ges=rho[0]*strengths[2]
    Gsn=rho[2]*strengths[7]
    L=1/Dtemporal(omega,alpha,beta)
    denominator=CharacteristicEquation(omega,k,alpha,beta,gamma,r,rho,strengths,t0)
    numerator=Ges*Gsn*(L**2)*np.exp(1j*omega*t0/2)
    h=numerator/denominator
    return h

def SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,strengths,t0,k0):
    p_a=0
    for m in range(16):
        for n in range(16):
            kx=2*np.pi*m/Lx
            ky=2*np.pi*n/Ly
            k=np.sqrt(kx**2+ky**2);
            p_a+=np.abs(TransferFunction(omega,k,alpha,beta,gamma,r,rho,strengths,t0))**2*np.exp(-k**2/k0)*2*np.pi/Lx*2*np.pi/Ly
    return p_a

def PoblationSpatialSpectrum(omega,G,phi,tau,alpha,beta,gamma=116,r=0.086):
    phi_a=0
    if len(G)==1:
        phi_a=G*phi*1/Dtemporal(omega,alpha,beta)*1/Dspatial(omega, gamma, r)*np.exp(1j*omega*tau)
    else:
        for m in range(len(G)):
            phi_a+=G[m]*phi[m,:]*1/Dtemporal(omega,alpha,beta)*1/Dspatial(omega, gamma, r)*np.exp(1j*omega*tau[m])
    return phi_a

def PoblationSpectrum(omega,G,phi,tau,alpha,beta):
    phi_a=0
    if len(G)==1:
        phi_a=G*phi*1/Dtemporal(omega,alpha,beta)*np.exp(1j*omega*tau)
    else:
        for m in range(len(G)):
            phi_a+=G[m]*phi[m,:]*1/Dtemporal(omega,alpha,beta)*np.exp(1j*omega*tau[m])
    return phi_a

def STDPspectrum(omega,ratesA,tp):
    
    
    A=np.zeros((len(ratesA),len(omega)),dtype=complex)
    B=np.zeros((len(ratesA),len(omega)),dtype=complex)
    for m in range(len(ratesA)):
        A[m,:]=ratesA[m]*tp*(1+1j*omega)/(1+(omega*tp)**2)
        B[m,:]=tp*(1-1j*omega)/(1+(omega*tp)**2)
    H=A+B
    return H
    

def integrand(phia,phib,Gamma,H):
    integ=np.real(phia*np.conjugate(phib)*np.conjugate(H)*np.conjugate(Gamma))
    return integ

def integral(integrand,omega,B=1):
    integrall=integrand[0]
    for i in range(1,len(integrand)):
        integrall+=integrand[i]*(omega[i]-omega[i-1])/(2*np.pi)
    return B*integrall

f=np.arange(0.01,50,0.05)
omega=f*(2*np.pi)
alpha=45
beta=186
# alpha=83
# beta=769
gamma=116
r=0.086
t0=0.085
k0=10
Lx=0.5  
Ly=0.5
#N2
# strengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3]    
#N3
# strengths=[6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.6e-3]

#SWS
strengths=[5.54e-3, -5.65e-3, 1.53e-3, 0.286e-3, 1.12e-3, 2.67e-3, -1.73e-3, 9.22e-3]

import matplotlib.pyplot as plt
phin=1
phin_std=4
k=0
rho=calculateRho(phin,strengths)
phi_e=TransferFunction(omega, k, alpha, beta, gamma, r, rho, strengths, t0)
Gee=strengths[0]*rho[0]
Gei=strengths[1]*rho[0]
Ges=strengths[2]*rho[0]
Gre=strengths[3]*rho[1]
Grs=strengths[4]*rho[1]
Gse=strengths[5]*rho[2]
Gsr=strengths[6]*rho[2]
Gsn=strengths[7]*rho[2]

phi_s0=PoblationSpectrum(omega,Gsn,np.ones_like(omega)*phin,0,alpha,beta)
phi_r0=PoblationSpectrum(omega,Grs,phi_s0,0, alpha,beta)
phi_i0=PoblationSpectrum(omega,np.array([Gee,Ges]),np.vstack((phi_e,phi_s0)),np.array([0,t0/2]), alpha, beta)
phi_e0=PoblationSpatialSpectrum(omega,np.array([Gee,Gei,Ges]),np.vstack((np.vstack((phi_e,phi_i0)),phi_s0)),np.array([0,0,t0/2]), alpha, beta)

phi_s0=phi_s0/np.abs(phi_s0)[0]
phi_r0=phi_r0/np.abs(phi_r0)[0]
phi_i0=phi_i0/np.abs(phi_i0)[0]
phi_e0=phi_e0/np.abs(phi_e0)[0]

phi_s1=PoblationSpectrum(omega,np.array([Gse,Gsr,Gsn]),np.vstack((np.vstack((phi_e0,phi_r0)),np.ones_like(omega)*phin)),np.array([t0/2,0,0]),alpha,beta)
phi_r1=PoblationSpectrum(omega,np.array([Gre,Grs]),np.vstack((phi_e0,phi_s1)),np.array([t0/2,0]),alpha,beta)
phi_i1=PoblationSpectrum(omega,np.array([Gee,Gei,Ges]),np.vstack((np.vstack((phi_e0,phi_i0)),phi_s1)),np.array([0,0,t0/2]), alpha, beta)
phi_e1=PoblationSpatialSpectrum(omega, np.array([Gee,Gei,Ges]), np.vstack((np.vstack((phi_e0,phi_i1)),phi_s1)),np.array([0,0,t0/2]), alpha, beta)

phi_s1=phi_s1/np.abs(phi_s1)[0]
phi_r1=phi_r1/np.abs(phi_r1)[0]
phi_i1=phi_i1/np.abs(phi_i1)[0]
phi_e1=phi_e1/np.abs(phi_e1)[0]


for m in range(8):
    phi_s2=PoblationSpectrum(omega,np.array([Gse,Gsr,Gsn]),np.vstack((np.vstack((phi_e1,phi_r1)),np.ones_like(omega)*phin)),np.array([t0/2,0,0]),alpha,beta)
    phi_r2=PoblationSpectrum(omega,np.array([Gre,Grs]),np.vstack((phi_e1,phi_s2)),np.array([t0/2,0]),alpha,beta)
    phi_i2=PoblationSpectrum(omega, np.array([Gee,Gei,Ges]), np.vstack((np.vstack((phi_e1,phi_i1)),phi_s2)),np.array([0,0,t0/2]), alpha, beta)
    phi_e2=PoblationSpatialSpectrum(omega, np.array([Gee,Gei,Ges]), np.vstack((np.vstack((phi_e1,phi_i2)),phi_s2)),np.array([0,0,t0/2]), alpha, beta)
    phi_e1=phi_e2/np.abs(phi_e2)[0]
    phi_i1=phi_i2/np.abs(phi_i2)[0]
    phi_s1=phi_s2/np.abs(phi_s2)[0]
    phi_r1=phi_r2/np.abs(phi_r2)[0]

#h=SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,strengths,t0,k0)

# plt.loglog(f,np.abs(phi_e)*f)
# plt.loglog(f,np.abs(phi_e2)*f)

ratesA=np.array([-1.01,-1.01,-1.25,-1,-0.4,-1.25,-1.25,-1.25])
tp=0.01


H=STDPspectrum(omega,ratesA,tp)
plt.figure()
plt.plot(np.abs(H[:,:].T))
plt.show()

Iee=integrand(phi_e,phi_e,1,H[0,:])
Iei=integrand(phi_e2,phi_i2,1,H[1,:])
Ies=integrand(phi_e2,phi_s2,1,H[2,:])


plt.figure()
#plt.plot(f,np.arcsinh(np.abs(phi_e2/np.sum(abs(phi_e2)))))
plt.plot(f,np.arcsinh(Iee))
#plt.plot(f,np.arcsinh(Iei))

plt.show()

print(integral(Iee,omega))
print(integral(Iei,omega))
print(integral(Ies,omega))
# plt.loglog(f,np.abs(phi_s2)*f)
# plt.loglog(f,np.abs(phi_r2)*f)
