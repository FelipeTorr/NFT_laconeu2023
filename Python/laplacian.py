import numpy as np
import matplotlib.pyplot as plt
import Data
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
import sys
from matplotlib import rc
rc('text', usetex=True)
def main():
	if len(sys.argv) > 1:
		#Load txt file
		filename = sys.argv[1]
		print('Loading and plotting '+filename+' ...')
		data, datazScore=Data.loadStoredData(filename)
		data=datazScore.to_numpy(dtype=float)
		nodes=256
		data=data[:,0:nodes]
		corr_e=np.dot(data.T,data)/np.var(data)
		#Laplacian = D-A
		A=np.zeros((nodes,nodes)) #binary adjacency matrix of more correlated nodes 
		D=np.zeros((nodes,nodes)) #Degree matrix
		k = int(sys.argv[2]) #Number of neighbors
		print(k)
		#A
		for u in range(nodes):
			neighbors=np.argsort(corr_e[u,:])
			A[u,neighbors[0:k]]=1
		
		for u in range(nodes):
			for v in range(nodes):
				if A[u,v]==1:
					A[v,u]=1
		#D
		for v in range(nodes):
			D[v,v]=np.sum(A[:,v])
		
		#Egigenvectors and eigenvalues
		eig_val,eig_vec=np.linalg.eig(D-A)
		norm = colors.Normalize(vmin=-0.5, vmax=0.5)
		normSum = colors.Normalize(vmin=0, vmax=2)
		fig1=plt.figure(constrained_layout=True,figsize=(9,9))
		gs1 = gridspec.GridSpec(4, 4, figure=fig1, width_ratios=[0.25, 0.25, 0.25, 0.25])
		for m in range(16):
			ax=fig1.add_subplot(gs1[m//4,m%4])		
			im=ax.imshow(np.reshape(eig_vec[:,m],(16,16)).T,cmap=cm.seismic,norm=norm)
			colorbar=plt.colorbar(im,ax=ax)

		plt.figure()
		plt.imshow(np.reshape(np.sum(np.abs(eig_vec[:,1:16]),axis=1),(16,16)).T,cmap=cm.hot_r,norm=normSum)
		plt.colorbar()
		plt.show()
	else:
		print ("Error: you must input the filename")
		sys.exit()
if __name__=='__main__':
	main()
