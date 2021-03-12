
class Conv_Operator:
	import numpy as np
	tol=10**-16
		
	def __init__(self,apply_elliptic_operator,order=2):
		self.order=order
		self.delta=lambda zeta : self.char_functions(zeta,order)
		self.apply_elliptic_operator=apply_elliptic_operator
	def get_integration_parameters(self,N,T):
		tol=self.tol
		dt=(T*1.0)/N
		L=N*2
		rho=tol**(1.0/(2*L))
		return L,dt,tol,rho

	def char_functions(self,zeta,order):
		if order==1:
			return 1-zeta
		else:
			if order==2:
				return 1.5-2.0*zeta+0.5*zeta**2
			else:
				if order ==3:
					#return 1-zeta**3
					return (1-zeta)+0.5*(1-zeta)**2+1.0/3.0*(1-zeta)**3
				else:
					print("Multistep order not availible")

	def get_zeta_vect(self,N,T):
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		import numpy as np
		Unit_Roots=np.exp(-1j*2*np.pi*(np.linspace(0,L,L+1)/(L+1)))
		Zeta_vect=self.delta( rho* Unit_Roots)/dt 
		return Zeta_vect

	def get_frequencies(self,N,T):
		import numpy as np
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		Unit_Roots=np.exp(-1j*2*np.pi*(np.linspace(0,L-1,L)/(L)))
		return rho*Unit_Roots

	def get_method_characteristics(self,method):
		import numpy as np
		import math
		if (method == "RadauIIA-2"):		
			c_RK=np.array([1.0/3,1])	
			A_RK=np.array([[5.0/12,-1.0/12],[3.0/4,1.0/4]])
			b_RK=np.array([[3.0/4,1.0/4]])	
		elif (method == "RadauIIA-3"):
			A_RK=np.array([[11.0/45-7*math.sqrt(6)/360, 37.0/225-169.0*math.sqrt(6)/1800 , -2.0/225+math.sqrt(6)/75],[37.0/225+169.0*math.sqrt(6)/1800,11.0/45+7*math.sqrt(6)/360,-2.0/225-math.sqrt(6)/75],[4.0/9-math.sqrt(6)/36,4.0/9+math.sqrt(6)/36,1.0/9]])
			c_RK=np.array([2.0/5-math.sqrt(6)/10,2.0/5+math.sqrt(6)/10,1])
			b_RK=np.array([4.0/9-math.sqrt(6)/36,4.0/9+math.sqrt(6)/36,1.0/9])

		m=len(A_RK[0,:])
		return A_RK,b_RK,c_RK,m	
	def format_rhs(self,rhs,m):
		try:
			N=(len(rhs[0,:]))/m
		except:
			N=(len(rhs))/m
			rhs_mat=np.zeros((1,N+1))
			rhs_mat[0,:]=rhs
			rhs=rhs_mat	
	
		return rhs,N

	def apply_RKconvol(self,rhs,T,show_progress=True,method="RadauIIA-2",cutoff=10**(-8)):
		import numpy as np
		import bempp.api
		## Step 1
		[A_RK,b_RK,c_RK,m]=self.get_method_characteristics(method)

		self.m=m	
		[rhs,N]=self.format_rhs(rhs,m)
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		rhs_fft=1j*np.zeros((len(rhs[:,0]),m*L))
		

		for stageInd in range(m):
			rhs_fft[:,stageInd:m*L:m]=self.scale_fft(rhs[:,stageInd:m*N:m],N,T)
		#Initialising important parameters for the later stage	
		s_vect=self.get_frequencies(N,T)

		dof=len(rhs[:,0])
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		Half=int(np.ceil(float(L)/2.0))

		## Step 2
		#phi_hat=1j*np.zeros((dof,L+1))
		normsRHS=np.ones(m*L)
		counter=0

		for j in range(0,m*L):
			normsRHS[j]=np.max(np.abs(rhs_fft[:,j]))
			if normsRHS[j]>cutoff:
				counter=counter+1

		#import matplotlib.pyplot as plt
		#plt.semilogy(normsRHS)
		#plt.show()
		print("CUTOFF :", cutoff)
		HalfL= int(np.ceil(float(L)/2.0))
		print("Amount of Systems needed: "+ str(counter))
		#Timing the elliptic systems
		import time
		start=0
		end=0

		rhsStages=1j*np.zeros((len(rhs_fft[:,0]),m))
		for j in range(0,HalfL+1):
			s=s_vect[j]
			deltaMatrix=np.linalg.inv(A_RK+s*1.0/(1-s)*np.ones((m,1))*b_RK)
			deltaEigs,T =np.linalg.eig(deltaMatrix)
		#	print("CONDITION T: ", np.linalg.cond(T))
			deltaEigs=deltaEigs/dt
			Tinv=np.linalg.inv(T)
#			print("NormT: ",np.linalg.norm(T))
			if show_progress:
				print("j:",j,"L:",L, "Time of previous iteration: " +str((end-start)/60), " MIN" )
			start=time.time()
			if j>0:

				relevantChange=False
				rhsStages=1j*rhsStages*0
				lhsStages=1j*lhsStages*0

				for stageInd in range(m):
					for sumInd in range(m):
						rhsStages[:,stageInd]=rhsStages[:,stageInd]+Tinv[stageInd,sumInd]*rhs_fft[:,m*j+sumInd]

					maxRHS=np.max(np.abs(rhsStages[:,stageInd]))
					#print("RHS MAX : ",maxRHS )
					if maxRHS>cutoff:
						relevantChange=True
						lhsStages[:,stageInd]=self.apply_elliptic_operator(deltaEigs[stageInd],rhsStages[:,stageInd])		
			else:
				relevantChange=True
				for sumInd in range(m):
					rhsStages[:,0]=rhsStages[:,0]+Tinv[0,sumInd]*rhs_fft[:,m*j+sumInd]
				first_eval=self.apply_elliptic_operator(deltaEigs[0],rhsStages[:,0])
				phi_hat=1j*np.zeros((len(first_eval),m*L))
				lhsStages=1j*np.zeros((len(phi_hat[:,0]),m))
				lhsStages[:,0]=first_eval
				for stageInd in range(1,m):
					for sumInd in range(m):
						rhsStages[:,stageInd]=rhsStages[:,stageInd]+Tinv[stageInd,sumInd]*rhs_fft[:,m*j+sumInd]
					lhsStages[:,stageInd]=self.apply_elliptic_operator(deltaEigs[stageInd],rhsStages[:,stageInd])		
			if relevantChange:	
				for stageInd in range(m):
					for sumInd in range(m):
						phi_hat[:,m*j+stageInd]=phi_hat[:,m*j+stageInd]+T[stageInd,sumInd]*lhsStages[:,sumInd]
			end=time.time()
#		for j in range(Half+1,L+1):
#			phi_hat[:,j]=np.conj(phi_hat[:,L+1-j])		
		
		## Step 3
		for freqInd in range(HalfL,L):	
			for stageInd in range(m):
				phi_hat[:,freqInd*m+stageInd]=np.conj(phi_hat[:,m*(L-freqInd)+stageInd])
		#print(phi_hat)		
		#print(phi_hat2)
		#print("ERR: ",np.abs(phi_hat-phi_hat2))
		#print("phi_hat ",phi_hat)
		#
		#print("phi_hat2 ",phi_hat2)
		#phi_sol=1j*phi_hat*0
		for stageInd in range(m):
			phi_hat[:,stageInd:m*L:m]=self.rescale_ifft(phi_hat[:,stageInd:m*L:m],N,T)
		#phi_sol[:,1:2*N+2:2]=self.rescale_ifft(phi_hat[:,1:2*N+2:2])
		phi_sol=phi_hat[:,:m*N]
		if len(phi_sol[:,0])==1:
			phi_sol=phi_sol[0,:]
		
		
		return phi_sol


	def apply_convol(self,rhs,T,show_progress=True,method="BDF2",cutoff=10**(-8)):
		import numpy as np
		import bempp.api
		if (method == "RadauIIA-2"):
			print(show_progress)
			return self.apply_RKconvol(rhs,T,show_progress=show_progress,method=method,cutoff=cutoff)

		## Step 1
		print("In apply_convol")
		try:
			N=len(rhs[0,:])-1
		except:
			N=len(rhs)-1	
			rhs_mat=np.zeros((1,N+1))
			rhs_mat[0,:]=rhs
			rhs=rhs_mat

		normsRHS=np.ones(N+1)
		for j in range(0,N+1):
			normsRHS[j]=np.max(np.abs(rhs[:,j]))

	
		#import matplotlib.pyplot as plt
		#plt.semilogy(normsRHS)
		#plt.show()
		
		rhs_fft=self.scale_fft(rhs)
		

#		import matplotlib.pyplot as plt
#		for j in range(0,len(rhs[:,0])):
#			print("rhs_fft : " , np.abs(rhs[j,:]))
#			plt.loglog(list(map(max,np.abs(rhs_fft[j,:]),10**-15*np.ones(len(rhs_fft[j,:])))))
#		plt.show()
#	

		#Initialising important parameters for the later stage	
		Zeta_vect=self.get_zeta_vect(N,T)
		print("ZETA0 :", Zeta_vect[0])
		dof=len(rhs[:,0])
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		Half=int(np.ceil(float(L)/2.0))

		## Step 2
		#phi_hat=1j*np.zeros((dof,L+1))
		normsRHS=np.ones(L+1)
		counter=0

		for j in range(0,Half+1):
			normsRHS[j]=np.max(np.abs(rhs_fft[:,j]))
			if normsRHS[j]>cutoff:
				counter=counter+1

		#import matplotlib.pyplot as plt
		#plt.semilogy(normsRHS)
		#plt.show()
		print("CUTOFF :", cutoff)
		
		print("Amount of Systems needed: "+ str(counter))
		#Timing the elliptic systems
		import time
		start=0
		end=0
		for j in range(0,Half+1):
			normsRHS[j]=np.max(np.abs(rhs_fft[:,j]))
			#print("normRHS:",normsRHS[j])
			if normsRHS[j]>cutoff:
				if show_progress:
					print("j:",j,"L:",str(Half), "Time of previous iteration: " +str((end-start)/60), " MIN" )
				start=time.time()
				if j>0:
					phi_hat[:,j]=self.apply_elliptic_operator(Zeta_vect[j],rhs_fft[:,j])
				else:
					first_eval=self.apply_elliptic_operator(Zeta_vect[j],rhs_fft[:,j])
					phi_hat=1j*np.zeros((len(first_eval),L+1))
					phi_hat[:,0]=first_eval
				end=time.time()
		for j in range(Half+1,L+1):
			phi_hat[:,j]=np.conj(phi_hat[:,L+1-j])		
		
		## Step 3
		import matplotlib.pyplot as plt

		phi_sol=self.rescale_ifft(phi_hat)
	
		phi_sol=phi_sol[:,:m*N]
		if len(phi_sol[:,0])==1:
			phi_sol=phi_sol[0,:]

		return phi_sol


	def scale_fft(self,A,N,T):
		import numpy as np
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		n_rows=len(A[:,0])
		n_columns=len(A[0,:])

		A=rho**(np.linspace(0,n_columns-1,n_columns))*A
		A_fft=np.fft.fft(np.concatenate((A,np.zeros((n_rows,L-n_columns))),axis=1))
		return A_fft

	def rescale_ifft(self,A,N,T):
		import numpy as np
		L,dt,tol,rho=self.get_integration_parameters(N,T)
		n_rows=len(A[:,0])
		n_columns=len(A[0,:])


		A_ift=np.real(np.fft.ifft(A))
		A_sol=np.real(rho**(-np.linspace(0,n_columns-1,n_columns))*A_ift[:,0:n_columns])
			
		return(A_sol)

