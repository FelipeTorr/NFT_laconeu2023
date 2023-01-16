class Population:
	#"""Population	class, to define paramaeters of	the	model and to tag results"""
	def __init__(self,name,N,L,v_matrix,q_matrix,phi_matrix,connections_matrix,Qmax=340,theta=0.01292,sigma_p=0.0038):
		#"""Constructor, requires	the	connections	array and the repective	input indexes
		#name: Population	name
		#Qinit: Initial firing rate
		#Qmax: maximum firing	rate
		#theta: firing rate threshold
		#sigma_p:	variance around	theta*sqrt(3)/pi"""
		self.name=name
		self.Qmax=Qmax
		self.theta=theta
		self.sigma_p=sigma_p
		self.connections=connections_matrix[(N-1)*L:N*L,(N-1)*L:N*L]
		self.V=v_matrix[(N-1)*L:N*L,:]
		self.Q=q_matrix[(N-1)*L:N*L,:]
		self.Phi=phi_matrix[(N-1)*L:N*L,:];

