#%% Baseline simulation of NeuralFields
import numpy as np
from solver import Solver
from model import Model

#%Common parameters for N2 and N3
	
#Spatial constants
gamma=116  #%(s^{-1}) wave damping factor
r=0.086  #%(m) axonal distance
t0=0.085 #%(s) cortico-thalamic delay
k=0
#%Frequency analysis parameters
f_final=50
omega=np.arange(0,2*np.pi*f_final+1,1)

#%%%% Simulation solution %%%%
#%%% No update connections: 0 % Update all connections: 1
plasticity_matrix=np.array([[0]]);
#%Population's characteristics
M=5
Nx=14 #%Lx paper=Nx*dx
Ny=14 #%Ly paper=Ny*dx
Lx=0.5 #%m
Ly=0.5 #%m
#%%% Activation function parameters (Sigmoid(x)=Qmax/(1+exp(-()))
Qmax=340
theta=0.01292
sigma_p=0.0038
#% Solver Time and Space parameters
h=1e-4 #%Step time
save_dt=1e-2 #%Sampling time
dx=Lx/Nx; #%0.017; #%30 nodes
final_Time=20 #%Simulation time
result_start_Time=3 #%time offset to avoid  
#%No stimulation
stimulation_Mode=0 #%0: no stimulation, 1: periodic stimulation, 2: in-phase stimulation, 3: random rate, 4:sinusoidal, 5:single pulse 
seed=0 #%Random number generator seed;
flag_update=0 #%If plasticity matrix is different from zero, this could avoid connections updating
	
#%% Sleep stage dependent parameters
strengths=np.zeros((4,8))
#%[vee, vei, ves, vre, vrs, vse, vsr, vsn]
#% N1
strengths[0,:]=np.array([5.73e-3, -6.39e-3, 0.24e-3, 0.87e-3, 0.52e-3, 6.07e-3, -1.46e-3, 28.52e-3])
#%N2 spindles
strengths[1,:]=np.array([3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3])
#%N3
strengths[2,:]=np.array([6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.58e-3])
#%Awake
strengths[3,:]=np.array([1.35e-3, -2.32e-3, 1.18e-3, 0.26e-3, 0.08e-3, 2.66e-3, -1.18e-3, 0.14e-3])

#%Rise and decay time
alpha=np.array([45,45,45,83],)
beta=np.array([186,186,186,769],)	

#%% Amplitudes of noise

asd=7.071067811865475e-05 #%Abeysuriya2014

if Nx*Ny>1:
	std_noise=np.sqrt(2*4*np.pi**3*asd**2/h/dx**2);
else:
	std_noise=np.sqrt(2.0*np.pi*asd**2/h)

print('Start \n');
for stage in range(1):		
#%%% Model instatiation and one dumb population (all population wit similar activation function) 
	
	modelEIRS=Model(strengths[stage:stage+1,:], plasticity_matrix, M, Nx, Ny,alpha[stage], beta[stage],gamma,r)
	modelEIRS.populations.append(modelEIRS.set_dumbPopulation(Qmax,theta,sigma_p))
	modelEIRS.set_t0(t0)

	#%%% Solver instatiation
	solver=Solver(modelEIRS,h,save_dt,dx,final_Time,result_start_Time,stimulation_Mode)
	solver.type_stim=0 #% Type auditory=0; tactile=1;
	solver.amplitude_stim=0
	solver.duration_stim=0 #%samples
	solver.frequency_stim=0 #%Hz
	#% Input noise
	solver.input_mean=0
	solver.input_std=std_noise
	solver.input_coherence=0
	filename='N'+str(stage)+'-baseline-python.csv'
	solver=solver.solve(seed,flag_update,filename);

print('The END \n');
