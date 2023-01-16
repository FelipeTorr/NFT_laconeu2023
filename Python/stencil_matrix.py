import numpy as np
def stencil_matrix(Nx,Ny,p2):
	W=Nx*Ny
	S=np.zeros((W,W))
	non_zero_indexes=np.zeros((W,4))
	for	m in range(Nx):#filas
		for n in range(Ny): #columnas
			#Left
			if n==0:
				non_zero_indexes[m*Ny+n,0]=m*Ny;
			else:
				non_zero_indexes[m*Ny+n,0]=m*Ny+n-1;	 
			

			#Up
			if m==0:
				non_zero_indexes[m*Ny+n,1]=(Nx-1)*Ny+n;	
			else:
				non_zero_indexes[m*Ny+n,1]=(m-1)*Ny+n;
			

			#Right
			if n==Ny-1:
				non_zero_indexes[m*Ny+n,2]=m*Ny;
			else:
				non_zero_indexes[m*Ny+n,2]=m*Ny+n+1;

			#Down
			if m==Nx-1:
				non_zero_indexes[m*Ny+n,3]=n;
			else:
				non_zero_indexes[m*Ny+n,3]=(m+1)*Ny+n;
					  
			
	
	
	for	j in range(W):
		for k in range(4):
			S[j,int(non_zero_indexes[j,k])]=p2
	return S
	