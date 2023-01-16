from population	import Population 
from model import Model	
from constants import * 
import time
import numpy as	np
from stencil_matrix import stencil_matrix

class Solver:
	#%% Constructor
	def __init__(self, model, h=1e-4,	dt_save=1e-2, dx=0.1, final_Time=100, result_start_Time=5, stimulation_Mode=0):
		if h>0:
			self.h=h
		else:
			print('Solver step time must be positive	non-zero\n')
			print('Set h=1e-4\n')
			self.h=1e4
		
		if dt_save>h:
			self.dt_save=dt_save
		else:
			print('Save time must be greather than solver step \n')
			print('Setting saving time to 100*h=	',100*h)
			self.dt_save=100*h
		
		
		if dx>0:
			self.dx=dx
		else:
			print('Delta space must be positive non-zero	\n')
			print('Set dx=0.01\n')
			self.dx=0.01
		
		
		self.second=np.ceil(1/h)
		self.samples_save=np.ceil(self.second*self.dt_save)
		self.model=model
		
		#Offset time to save results must	be zero	or greater than	one
		if result_start_Time==0:
			self.result_start_time=0
		elif result_start_Time>=1:
			self.result_start_time=result_start_Time
		else:
			self.result_start_time=1
			print('Start time offset	must be	superior than one second or	zero, set =1\n')
		
		
		if final_Time>0:
			self.final_time=final_Time
		else:
			print('Final time must be positive non-zero \n')
			print('Final time set to	result_start_time+1: ',result_start_Time+1)
			self.final_time=self.result_start_time+1

		
		#Stimulation mode, could be
		#0: no stimulation, 1:constant rate pulse, 2:phase-locked
		#3:random	rate pulse,	4:sinusoidal, 5: single	pulse
		self.stimulation_mode=stimulation_Mode
		#pulse duration
		self.duration_stim=1
		#pulse amplitude
		self.amplitude_stim=1e-7
		#constant	rate frequency or mean of random rate
		self.frequency_stim=0.5
		#To change between pulses
		self.type_stim=0
		#input noise noise std
		self.input_std=1e-5
		#input noise mean
		self.input_mean=0
		#Noise for each node (0) or for all nodes	(1)
		self.input_coherence=0
		#random rate std
		self.frequency_std_stim=0
		#mean	V_e	voltage
		self.mean_v=0
		#mean	Q_e	firirng	response
		self.mean_q=0
		#	Threshold for slow oscillations	detection
		self.SO_th=-40e-6
		#	Flag to	save mean_v	and	mean_q values
		self.save_means_flag=0
		#	Phase targe	for	pahase-locked stimulation
		self.target_phase_stim=0
		#	Distribution of	random inputs 0: uniform, 1:normal
		self.random_type=1
		
		self.dt2on12=h**2/12
		self.dfact=self.dt2on12*self.model.gamma**2
		self.dt2ondx2=h**2/dx**2;
		self.p2=self.dt2ondx2*self.model.range_e**2*self.model.gamma**2
		self.tenminus4p2=10-4*self.p2
		self.twominus4p2=2-4*self.p2
		self.expfactneg=np.exp(-h*self.model.gamma)
		self.expfactpos=np.exp(h*self.model.gamma)
		self.Sp2=stencil_matrix(self.model.Nx,self.model.Ny,self.p2)
		self.I=np.eye(self.model.L)
	
		
	#%% Auxiliar functions
	#Firing funtions
	def sigmoid(self,V,Qmax,theta,sigma_p):
		#"""Sigmoid function for firing response"""
		firing_response=Qmax/(1+np.exp(-(V-theta)/sigma_p));
		return firing_response
	
	def linear(self,V,Qmax,Q0,a):
		#"""Sigmoid function for firing response"""
		firing_response=Q0+a*V
		if firing_response>Qmax:
			firing_response=Qmax
		return firing_response
	
	
	#Shift	samples
	#Circshift	is better in MATLAB
	def shift_sample(a):
		size_a=np.shape(a)
		if size[0]>size[0]:
			a[1:size_a[0],0]=a[1:size_a[0]-1,0]
		else:
			a[0,1:size_a[1]]=a[0,0:size_a[1]-1]
			
	#%%%% Solve methods %%%%%%%%%%%
	#%% Principal encapsulating method
	def solve(self, seed,update_flag,filename):
		print('Starting solver...\n');
		fileID = open(filename,"w")
		print('Writing results in ', filename);
		#%%% Additional constants
		
		delay_samples=np.floor(self.model.t0/(2*self.h))
		connections_update_flag=update_flag #%%TODO: update of connection	weights
		total_time_samples=np.ceil(self.final_time/self.h)
		j_connection=1 #connection matrix	updating index or seconds counter
		marker=0
		if self.frequency_stim>0:
			stimulus_time=np.ceil(1/(self.h*self.frequency_stim));
		else:
			stimulus_time=0
		pre_stimulus=np.ceil(2/self.h);
		stim_flag=0;
		end_stim=0;
		#%Output vectors and matrices
		#%phi_matrix=zeros(self.model.W,delay_samples+1);
		v_matrix=np.zeros((self.model.W,1));
		s1,s2=np.shape(self.model.connection_matrix)
		history_connections=np.zeros((s1,s2,int(np.ceil(self.final_time))));
		[q_pre,phi_pre,vab_matrix_prev]=self.presolver(seed);
		print('Pre-solver stage	completed\n')
		#Change in the time order, first element now corresponds to
		#last	time index
		solverTime=time.time()
		phi_matrix=np.copy(phi_pre)
		q_matrix=np.copy(q_pre)
		 
		solver_loops=np.ceil(total_time_samples/delay_samples);
		print('Start of saving data simulation\n')
		for loop in range(int(solver_loops)):
		#Random input
			if self.input_coherence==1: #All	nodes same input
				np.random.seed(seed+loop)
				rnd_input=np.random.randn(1,int(delay_samples)+1)*self.input_std+self.input_mean
				for	j in range(self.model.W-self.model.L,self.model.W):
					phi_matrix[j,:]=np.copy(rnd_input)
					#voltage at input are zero	or calculated as S^{-1}(Q_n)
						 
			else: #Different	input for each node
				for	j in range(self.model.W-self.model.L,self.model.W):
					np.random.seed(seed+self.model.L*j+loop)
					rnd_input=np.random.randn(1,int(delay_samples)+1)*self.input_std+self.input_mean
					phi_matrix[j,:]=np.copy(rnd_input)
					#voltage at input are zero	or calculated as S^{-1}(Q_n)
					
			for	n_loop in range(int(delay_samples)+1):
				#******************************************
				v_matrix,q_matrix,phi_matrix,vab_matrix_prev=self.step(v_matrix,q_matrix,phi_matrix,vab_matrix_prev)
				#******************************************

				#***********Stimulation start **********************
				current_time=(loop-1)*(delay_samples+1)+n_loop	# in samples
				#%%% %Only	White noise	%%%
				if	self.stimulation_mode==0:
					marker=0
					delta=0
				#%%% Pulse	constant rate %%%
				elif self.stimulation_mode==1:
					#%start of stimulation if	there is enough	time and another
					#%stimulus is	not	present
					if (current_time-pre_stimulus)%stimulus_time==0 and current_time+self.duration_stim<total_time_samples and stim_flag==0:
						stim_flag=1
						end_stim=current_time+self.duration_stim
						marker=1

					#	Add	stimulus amplitude to the input	population activity
					if stim_flag==1 and current_time<=end_stim:
						phi_matrix[self.model.W-self.model.L:self.model.W,0]=phi_matrix[self.model.W-self.model.L:self.model.W,0]+self.amplitude_stim;
					if current_time==end_stim:
						stim_flag=0;
						delta=0;

				#%%% Pulse	random rate	%%%
				elif self.stimulation_mode==2:
					if (current_time-pre_stimulus)% stimulus_time==0 and current_time+self.duration_stim<total_time_samples	and stim_flag==0:
						stim_flag=1
						end_stim=current_time+self.duration_stim
						marker=1
						#Change next	time
						if self.random_type==1: #%uniform
							stimulus_time=np.ceil(1/(self.h*(np.random.rand(1)*self.frequency_stim)))#Uniform betwwen 0	and	frequency_stim
						elif self.random_type==2:	#%normal
							stimulus_time=np.ceil(1/(self.h*np.max([self.frequency_std_stim*np.random.rand(1,1)+self.frequency_stim,0]))) #%normal and relu

					#Add stimulus	amplitude to the input population activity
					if stim_flag==1 and current_time<=end_stim:
						phi_matrix[self.model.W-self.model.L:self.model.W,0]=phi_matrix[self.model.W-self.model.L:self.model.W,0]+self.amplitude_stim;
						if current_time==end_stim:
							stim_flag=0


				#%%% Phase	locked stimulation %%%	
				elif self.stimulation_mode==3:
				#Stim	marker
					stim=0			   
					delta[current_time]=0 #while	the	first second of	stimulation
					#Phase detection	algorithm
					if	current_time>self.second:
						previous_delta[1::]=previous_delta[0:-1];
						previous_delta[0]=np.mean(v_matrix[0:self.L,0])-self.mean_v
						#filtering
						previous_delta_filter=core.filtrar(b_delta,a_delta,previous_delta,previous_delta_filter)	#%filtrar

						origen_eeg[1:3]=origen_eeg[0:2]
						origen_eeg[1]=filtro_delta
						shiftphase_eeg[1]=shiftphase_eeg[0]
						shiftphase_eeg[0]=(origen_eeg[2]-origen_eeg[0])
						x_now[1:3]=x_now[0:2]

						#Rectifiying
						if filtro_delta>0:
							x_now[0]=filtro_delta;
						else:
							x_now[0]=-filtro_delta;

						#Filtro pasabajo 0.5Hz Buttherworth 4	orden
						envolvente,previous_envelope_filter=core.filtrar(b_envelope,a_envelope,x_now,previous_envelope_filter)
						if envolvente<ENVELOPE_OFFSET: #%Envolvente menor	a 10uV
							envolvente=ENVELOPE_OFFSET #%Valor mínimo 10	uV
						else:
							envolvente+=ENVELOPE_OFFSET #%sumar	el margen pues la envolvente
							#no necesariamente es el maximo

						#delta(i)=filtro_delta;
						#Normalización
						eeg_normalizado=origen_eeg[0]/envolvente;
						eeg_normalizado_90=shiftphase_eeg[0]/envolvente;
						#Calculo de la fase a	partir del arco	tangente
						fase=atan(-eeg_normalizado_90/(eeg_normalizado+1e-12));

						#Para	pasar del dominio [-pi,	pi]	al dominio [0, 2pi]
						if origen_eeg[0]>0 and origen_eeg[1]<0: #seno	subida
							#aumento_fase=pi/2; python
							aumento_fase=pi/2;
						elif origen_eeg[1]>0 and origen_eeg[0]<0:	#seno bajada
							#aumento_fase=3/2*pi; python
							aumento_fase=3*pi/2;
						elif shiftphase_eeg[0]>0 and shiftphase_eeg[1]<0:	#coseno	subida
							#aumento_fase=5/2*pi; python
							aumento_fase=3/2*pi;
						elif shiftphase_eeg[1]>0 and shiftphase_eeg[0]<0:	#coseno	bajada
							#aumento_fase=pi/2; python
							aumento_fase=pi/2;

						#sumar el	desfase	correspondiente	segun el
						#cruce por cero
						fase=fase+aumento_fase;
						delta=filtro_delta;
						#Si la fase resulta mayor	a 2pi mantener en esa cota
						if fase>2*pi:
							phase=2*pi
						elif fase<0:
							phase=0
						else:
							phase=fase


						if (previous_delta(1)<self.SO_th	and flag_hold==0):
							flag_search=1						
						if (flag_search==1 and phase>=self.target_phase_stim-EPSILON_ANGLE and phase<=self.target_phase_stim+EPSILON_ANGLE and filtro_delta>0):
							if flag_shot==0:
								flag_shot=1    #Perform a	pulse of the stimulus
								flag_hold=1    #A stimulus consists	in two pulses, this	flag guarantee the delivery	of both
								flag_search=0  #Deactivate for not search the angle	again until	the	first pulse	is delivered

						if flag_shot==1:
							flag_shot=0
							if STIM_SHAM==1:
								#If	STIM radio button is selected
								if flag_pause==1:
									#Similar to SHAM condition
									if	pulses==0:
										pulses=1;
										stim=2
									else:
										pulses=3;
										stim=4;

								else:
									if	STIM_TYPE==0:
										print("Auditory stimulus \n")
										phi_matrix[W-L:W,0:int(np.ceil(self.duration_stim))]=phi_matrix[W-L:W,0:int(np.ceil(self.duration_stim))]+self.amplitude_stim	
									else:
										#TODO:	tactile	stimulation
										#Hint:	could be to	play a different sound to activate the actuator

										print("Tactile stimulus	\n")

									if pulses==0: #Control of markers and delivery of two pulses
										pulses=1
										stim=1
									else:
										pulses=3
										stim=3
										flag_pause=1
										time_pause=current_time;
										COUNT_STIMULI+=1
										#%print("start pause \n")

							else:
								#SHAM
								if pulses==0:
									pulses=1
									stim=2
								else:
									pulses=3
									stim=4

						if (flag_hold==1 and phase>pi):
							if pulses==1:
								if (filtro_delta<self.SO_th):
									flag_search=1
									pulses=2

							elif pulses==3:
								flag_hold=0
								pulses=0

						if (current_time-time_pause)>PAUSE_TIME*self.second and	flag_pause:
							flag_pause=0
							flag_hold=0
							pulses=0
							#%print("end	pause \n")

						marker=stim;

				#%%%	Sinusoidal input %%%	
				elif	self.stimulation_mode==4:
					#%start of stimulation if	there is enough	time and another
					phi_matrix[self.model.W-self.model.L:self.model.W,0]=phi_matrix[self.model.W-self.model.L:self.model.W,0]+self.amplitude_stim*np.sin(2*np.pi*self.frequency_stim*self.h*(current_time+1));
								  
				
				#%%%%%%%%%% Stimulation end %%%%%%%%%%%%%%%
				#%%%	Connections	update
				if (current_time%self.second)==0:
					if connections_update_flag==1:
						self.model.connection_matrix=self.update_conections(self.model.connection_matrix)
						history_connections[:,:,j_connection]=self.model.connection_matrix[:,:]
					
					print('Simulated',str(j_connection), ' s of ', filename)
					print('Phi: ',phi_matrix[[1,4*self.model.L+1],1])
					j_connection+=1 #one second more
					
				
				#%%%% Save data	if dt_save is accomplished
				if (current_time%self.samples_save)==0:
					#%Excitatory nodes expected	in (1:L) positions
					for m_save in range(self.model.L):
						fileID.write('%.12e\t' % phi_matrix[m_save,0])
						fileID.write('%.12e\t' % v_matrix[m_save,0])
					if self.stimulation_mode==3:
						fileID.write('%.12e\t\t' % delta)
						fileID.write('%.12e\t \n' % marker)
					else:
						fileID.write('%.12e\t \n' % marker)
					marker=0;
				
			  #end for delay_samples
		#end for solver_loops
		print('Solver complete:	');
		print('Elapsed time', time.time()-solverTime)
		fileID.close()
		
	#end solver function
	
	#%% Set values
	def set_meanValues(self,mean_v,mean_q,SO_th):
		self.mean_v=mean_v;
		self.mean_q=mean_q;
		self.SO_th=SO_th;
	
	
	#%% Pre-solver%%%
	def presolver(self,seed):
		time_presolver=np.ceil(self.result_start_time);
		W=self.model.W
		L=self.model.L
		A=self.model.A
		B=self.model.B
		connection_matrix_delay=self.model.connection_matrix_delay
		connection_matrix_nodelay=self.model.connection_matrix_nodelay
		Qmax=self.model.populations[0].Qmax
		theta=self.model.populations[0].theta
		sigma_p=self.model.populations[0].sigma_p
		delay_samples=np.floor(self.model.t0/(2*self.h))
		vab_matrix=np.zeros((self.model.W-self.model.L,2))
		meanv=0
		meanq=0
		#%%% Execute the solver for time_presolver seconds before	save results
		
		#Firing response vector
		q_matrix=np.zeros((int(self.model.W),3))
		#Potentials vector
		phi_matrix=np.zeros((int(self.model.W),int(delay_samples)+1))
		
		if self.input_coherence==1: #All nodes same input
			np.random.seed(seed)
			rnd_input=np.random.randn(1,self.second*time_presolver)*self.input_std+self.input_mean
			for	j_n in range(self.model.W-self.model.L,self.model.W):
				phi_matrix[j_n,1]=rnd_input[1,1]
				
		else: #Different input for each	node
			np.random.seed(seed)
			rnd_input=np.random.randn(int(self.model.L),int(self.second*time_presolver))*self.input_std+self.input_mean
			for	j_n	in range(self.model.L):
				phi_matrix[j_n+self.model.W-self.model.L,1]=rnd_input[j_n,1]
				#voltage	at input are zero or calculated	as S^{-1}(Q_n)
				
		
		#Matrix calculations
		print('Presolver: ')
		new1=time.time()
		for i_prev in range(int(np.ceil(self.second*time_presolver))):
			#Soma voltages	[V_a] at time i
			v_aux=connection_matrix_nodelay*phi_matrix[:,0]+connection_matrix_delay*phi_matrix[:,int(delay_samples)];
			#Soma voltages	at time	i+1	((W-L+1:end) terms are zero)
			vab_matrix=vab_matrix+self.h*(np.matmul(vab_matrix,A)+np.matmul(v_aux[0:W-L,0:1],B));
			
			#Shift elements (time consuming operation)
			q_matrix=np.roll(q_matrix,1,axis=1);
			phi_matrix=np.roll(phi_matrix,1,axis=1);
			#Firing responses [Q_a] for next	time i+1
			q_matrix[0:W-L,0]=self.sigmoid(vab_matrix[:,0],Qmax,theta,sigma_p);
			
			#Spatial	propagation: Excitatory	population nodes expected in (1:L)
			phi_matrix[0:L,1]=self.propagate(q_matrix[0:L,:],phi_matrix[0:L,1:3])
			#Another	populations	maps the firing	response
			phi_matrix[L:W-L,0]=q_matrix[L:W-L,0]
			
			minv=np.min(vab_matrix[0:L,0])
			#Noise input
			if self.input_coherence==1:
				for	j_n	in range(self.model.L):
					phi_matrix[self.model.W-self.model.L+j_n,0]=rnd_input[0,i_prev];
				
			else:
				for j_n in	range(self.model.L):
					phi_matrix[self.model.W-self.model.L+j_n,0]=rnd_input[j_n,i_prev];
			
			if	self.save_means_flag==1:
				#Means of coritcal population
				meanv=meanv+np.mean(vab_matrix[0:L,0])/(self.second*time_presolver-1);
				meanphi=meanphi+np.mean(phi_matrix[0:L,0])/(self.second*time_presolver-1);
				minv_candidate=np.min(vab_matrix[0:L,0]);
				if minv_candidate<minv:
					minv=minv_candidate
				
				SOth=0.75*(minv-meanv)
				self.set_meanValues(self,meanv,meanq,SOth)
				
		#end for
		print('elapsed time: ',time.time()-new1)
		return q_matrix,phi_matrix,vab_matrix
	#Presolver	function end
	
	#%% Step of the solver
	def step(self,v_matrix,q_matrix,phi_matrix,vab_matrix_prev):
		#Model characterisitics
		W=self.model.W
		L=self.model.L
		A=self.model.A
		B=self.model.B
		connection_matrix_delay=self.model.connection_matrix_delay
		connection_matrix_nodelay=self.model.connection_matrix_nodelay
		Qmax=self.model.populations[0].Qmax
		theta=self.model.populations[0].theta
		sigma_p=self.model.populations[0].sigma_p
		#Solver
		#Soma	voltages [V_a] at time i
		#	index 1: last time value, index	end: time-t0/2 sample
		#tic()
		v_aux=connection_matrix_nodelay*phi_matrix[:,0]+connection_matrix_delay*phi_matrix[:,-1]
		#fprintf('step 1:	')
		#toc()
		#Soma	voltages at	time i+1 ((W-L+1:end) terms	are	zero)
		#tic()
		vab_matrix_post=vab_matrix_prev+self.h*(np.dot(vab_matrix_prev,A)+np.dot(v_aux[0:W-L,0:1],B))
		#fprintf('step 2:	')
		#toc()
		#	Shift in time
		#tic()
		q_matrix=np.roll(q_matrix,1,axis=1)
		phi_matrix=np.roll(phi_matrix,1,axis=1)
		v_matrix[0:W-L,0]=vab_matrix_post[:,0] #%Store soma voltages (squeeze(vab))
		#fprintf('step 3:	')
		#toc()
		#%Firing responses [Q_a] for next	time i+1
		#tic()
		q_matrix[0:W-L,0]=self.sigmoid(v_matrix[0:W-L,0],Qmax,theta,sigma_p)
		#fprintf('step 4:	')
		#toc()
		#Spatial propagation
		#tic()
		phi_matrix[0:L,0]=self.propagate(q_matrix[0:L,0:3],phi_matrix[0:L,1:3])
		phi_matrix[L:W-L,0]=q_matrix[L:W-L,0];
		#fprintf('step 5:	')
		#toc()
		return v_matrix,q_matrix,phi_matrix,vab_matrix_post
	
	
	
	#%% Spatial propagation
	def propagate(self,q_matrix_3,phi_matrix_2):
		#q_matrix_3 last 3 values
		#phi_matrix_3	last 2 values
		#Note	the	different time indexes
		#index 3:	n-1, index 2: n, index 1: n+1
		#gamma \gamma, damping factor
		#range r_e, axonal range
		#deltat, deltax spatiotemporal discretization
		#Lx horizontal nodes
		#Ly vertical nodes
		#phi_matrix=zeros(Lx*Ly,1);
		#if deltat>gamma/2 &&	deltax>range/2
		drive=self.dfact*(np.dot((np.dot(self.tenminus4p2,self.I)+self.Sp2),q_matrix_3[:,1])+self.expfactpos*np.dot(self.I,q_matrix_3[:,0])+self.expfactneg*np.dot(self.I,q_matrix_3[:,2]))
		phi_matrix=self.expfactneg*(np.dot((np.dot(self.twominus4p2,self.I)+self.Sp2),phi_matrix_2[:,0])-self.expfactneg*np.dot(self.I,phi_matrix_2[:,1])+drive)
		#else:
		#fprintf('gamma and range must	be at maximum the double of	delta t	and	delta x, respectively.')
		return phi_matrix
	

	
	#%% Update	connections	strengths
	def update_connections(self, connection_matrix):
		new_connection_matrix=connection_matrix
		new_connection_matrix[connection_matrix>0]=connection_matrix[connection_matrix>0]+1e-9*self.model.plasticity_matrix[connection_matrix>0]
		new_connection_matrix[connection_matrix<0]=connection_matrix[connection_matrix<0]-1e-9*self.model.plasticity_matrix[connection_matrix<0]
		return new_connection_matrix
