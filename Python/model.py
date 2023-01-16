import numpy as np
from population	import Population 

class Model:
	def __init__(self,initial_weights, plasticity_matrix, M=5, Nx=1, Ny=1,alpha=45, beta=186,gamma=116,range_e=0.086):	  
		if M>=2: #At least two populations
			self.M=M
			print('Model: Set populations quantity to ',str(M))
		else:
			self.M=5
			print('Model: Set populations quantity to 2 \n')
		
		if (Nx>=1	and	Nx<50):	#At	least one node per row and no more than	50 nodes per row
			self.Nx=Nx
			print('Model: Set nodes per row to ',str(Nx))
		else:
			self.Nx=1
			print('Model: Set nodes per row to 1	\n')
		
		if (Ny>=1	and	Ny<50):	#At	least one node per column and no more than 50 nodes	per	column
			self.Ny=Ny
		else:
			self.Ny=1
			print('Model: Set nodes per column to 1 \n')
		
		#Calculate W
		self.L=self.Nx*self.Ny
		self.W=self.M*self.L
		#Delay matrix
		delay_aux_matrix=np.zeros((5,5));
		delay_aux_matrix[1,4]=1;
		delay_aux_matrix[2,4]=1;
		delay_aux_matrix[4,1]=1;
		delay_aux_matrix[3,1]=1;
		self.delay_matrix=self.expand_connections_strengths(delay_aux_matrix);
		
		size_ic=np.shape(initial_weights)
		if size_ic[0]==self.W and size_ic(1)==self.W:
			#If the input is	a convenient square	matrix
			self.connection_matrix=initial_weights;
			
			print('Model: Copied all	connections	values	\n');
		elif (size_ic[0]==1 and size_ic[1]==8):
			#By default use of EIRS model
			eirs_connections=self.eirs_connections_strengths(initial_weights);
			
			self.connection_matrix=self.expand_connections_strengths(eirs_connections);
			self.connection_matrix_delay=self.connection_matrix*self.delay_matrix;
			self.connection_matrix_nodelay=self.connection_matrix*(1-self.delay_matrix);
			print('Model: EIRS model \n');
		else:
			
			#A square matrix of the model connections is	expected (same value for all nodes of one population)
			self.connection_matrix=self.expand_connections_strengths(initial_weights);
			print('Model: Expanded connection wights to a matrix \n');
		
		
		size_plast=np.shape(plasticity_matrix)
		if (size_plast[0]==self.W and size_plast[1]==self.W):
			#If the input is	a convenient square	matrix
			self.plasticity_matrix=plasticity_matrix
			print('Model: Update of selected connections \n');
		elif plasticity_matrix==1:
			#If the input is	1, all connections are plastic
			self.plasticity_matrix=np.ones((self.W,self.W))
			print('Model: Update of all connections \n')
		else:
			#If the input is	0, (any	another	in fact), all connections
			#are	fixed.
			self.plasticity_matrix=np.zeros((self.W,self.W))
			print('Model: No connections update \n')
	
		#4. Differential equations of	synaptodendritic potentials
		self.A=np.transpose(np.array([[0,	1],	[-alpha*beta, -alpha-beta]]))
		self.B=np.transpose(np.array([[0],[alpha*beta]]));
		self.populations=[]
		self.t0=0.085	#seconds
		
		
		self.gamma=gamma;
		self.range_e=range_e;
		print('Model: Instantiatiation completed \n')
		
		
	def expand_connections_strengths(self, strengths_matrix):
		#"""Expands each strenght	between	populations	as an intra-population identity	matrix """
		identity=np.eye(self.L);
		connections_matrix=np.zeros((self.W,self.W))
		size_M=np.shape(strengths_matrix)
		for i in range(size_M[0]):
			for j in range(size_M[1]):
				connections_matrix[i*self.L:(i+1)*self.L,j*self.L:(j+1)*self.L]=identity*strengths_matrix[i,j];
		output_matrix=connections_matrix;
		return output_matrix
		
		
	def eirs_connections_strengths(self, eirs_strengths):
		#"""Construct	a matrix with the relevant strenghts"""
		#%%% EIRS	weights	(in	some papers	the	gains are indicated	Gab=rho*vab, and rho is	the	slope of the sigmoid ~340)
		#	[vee, vei, ves,	vie, vii, vis, vir,	vre, vrs, vse, vsr,	vsn]
		#	vee==vie
		#	vei==vii
		#	ves==vis
		#	[vee, vei, ves,	vre, vrs, vse, vsr,	vsn]
		strengths_matrix=np.zeros((5,5))
		strengths_matrix[0,0]=eirs_strengths[0,0]
		strengths_matrix[0,1]=eirs_strengths[0,1]
		strengths_matrix[0,3]=eirs_strengths[0,2]
		strengths_matrix[1,0]=eirs_strengths[0,0]
		strengths_matrix[1,1]=eirs_strengths[0,0]
		strengths_matrix[1,3]=eirs_strengths[0,2]
		strengths_matrix[2,0]=eirs_strengths[0,3]
		strengths_matrix[2,3]=eirs_strengths[0,4]
		strengths_matrix[3,0]=eirs_strengths[0,5]
		strengths_matrix[3,2]=eirs_strengths[0,6]
		strengths_matrix[3,4]=eirs_strengths[0,7]
			
		return strengths_matrix
	
	def set_dumbPopulation(self,Qmax,theta,sigma_p):
		#One dumb	population in just necessary now to	apply the matrix
		#calculation way
		newPopulation=Population("dumb",1,1,np.array([[0]]),np.array([[0]]),np.array([[0]]),self.connection_matrix,Qmax,theta,sigma_p)
		print('Model: Dumb population created \n')
		return newPopulation
	
	
	def setPopulations(self,Qmax,theta,sigma_p,v_matrix,q_matrix,phi_matrix):
		#Set the outputs to a	populations	array to understand	them
		#better
		#Assuming	EIRS
		for i, name in zip(range(self.M),['E','I','R','S','N']):
			newPopulation=Population(name,i,self.L,v_matrix,q_matrix,phi_matrix,self.connection_matrix,Qmax,theta,sigma_p)
			self.populations.append(newPopulation)
	def set_t0(self,t0):
		self.t0=t0
