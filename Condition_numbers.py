import numpy as np
import bempp.api
import math
from RKconv_op import *

def create_timepoints(c,N,T):
	m=len(c)
	time_points=np.zeros((1,m*N))
	for j in range(m):
		time_points[0,j:m*N:m]=c[j]*1.0/N*np.ones((1,N))+np.linspace(0,1-1.0/N,N)
	return T*time_points

def create_rhs(grid,N,T,m):
	#grid=bempp.api.shapes.sphere(h=dx)
	if (m==2):
		c_RK=np.array([1.0/3,1])
	if (m==3):
		c_RK=np.array([2.0/5-math.sqrt(6)/10,2.0/5+math.sqrt(6)/10,1])

	from bempp.api.operators.boundary import maxwell
	
	RT_space = bempp.api.function_space(grid,"RT",0)
#	curl_space=bempp.api.function_space(grid,"RBC",0)
	#from bempp.api.operators.boundary.sparse import identity as ident
#	id1 = ident(div_space,div_space,curl_space).weak_form()
#	print("CONDITION NUMBER : ", np.linalg.cond(bempp.api.as_matrix(id1).todense()))
	dof=RT_space.global_dof_count
	print(" DOF: ", dof)
	rhs=np.zeros((dof+dof,N*m))
	curls=np.zeros((dof,N*m))
	time_points=create_timepoints(c_RK,N,T)
	for j in range(m*N):
		t=time_points[0,j]
		def incident_field(x):
			return np.array([np.exp(-50*(x[2]-t+1)**2), 0. * x[2], 0. * x[2]])
			#return np.array([np.exp(-200*(x[2]-t+2)**2), 0. * x[2], 0. * x[2]])

		def tangential_trace(x, n, domain_index, result):
			result[:] = np.cross(n,np.cross(incident_field(x), n))

		def curl_trace(x,n,domain_index,result):
			curlU=np.array([ 0. * x[2],-100*(x[2]-t+1)*np.exp(-50*(x[2]-t+1)**2), 0. * x[2]])
			result[:] = np.cross(curlU , n)

		curl_fun = bempp.api.GridFunction(RT_space, fun=curl_trace,dual_space=RT_space)
		trace_fun= bempp.api.GridFunction(RT_space, fun=tangential_trace,dual_space=RT_space)
		rhs[0:dof,j]=trace_fun.coefficients	
		curlCoeffs=curl_fun.coefficients
		if np.linalg.norm(curlCoeffs)>10**-9:
			curls[0:dof,j]=curlCoeffs

		#print("RHS NORM :", np.linalg.norm(trace_fun.coefficients))

	def sinv(s,b):
		return s**(-1)*b
	IntegralOperator=Conv_Operator(sinv)
	def HarmonicImpedance(s,b):
		return 0.1*s**(0.5)*b
	TimeImpedance=Conv_Operator(HarmonicImpedance)	

	if (m==2):
		curls=IntegralOperator.apply_RKconvol(curls,T,method="RadauIIA-2",show_progress=False)
		ZptNeuTrace=TimeImpedance.apply_RKconvol(curls,T,method="RadauIIA-2",show_progress=False)
	if (m==3):
		curls=IntegralOperator.apply_RKconvol(curls,T,method="RadauIIA-3",show_progress=False)
		ZptNeuTrace=TimeImpedance.apply_RKconvol(curls,T,method="RadauIIA-3",show_progress=False)


	rhs[0:dof,:]=np.real(ZptNeuTrace)-rhs[0:dof,:]
	return rhs

def harmonic_calderon(s,b,grid,points):
	OrderQF = 8

	#tol= np.finfo(float).eps
	bempp.api.global_parameters.quadrature.near.max_rel_dist = 2
	bempp.api.global_parameters.quadrature.near.single_order =OrderQF-1
	bempp.api.global_parameters.quadrature.near.double_order = OrderQF-1
	
	bempp.api.global_parameters.quadrature.medium.max_rel_dist =4
	bempp.api.global_parameters.quadrature.medium.single_order =OrderQF-2
	bempp.api.global_parameters.quadrature.medium.double_order =OrderQF-2


	bempp.api.global_parameters.quadrature.far.single_order =OrderQF-3
	bempp.api.global_parameters.quadrature.far.double_order =OrderQF-3

	bempp.api.global_parameters.quadrature.double_singular = OrderQF
	bempp.api.global_parameters.hmat.eps=10**-4

	bempp.api.global_parameters.hmat.admissibility='strong'
###    Define Spaces
	NC_space=bempp.api.function_space(grid, "NC",0)
	RT_space=bempp.api.function_space(grid, "RT",0)
		
	elec = -bempp.api.operators.boundary.maxwell.electric_field(RT_space, RT_space, NC_space,1j*s)


	magn = -bempp.api.operators.boundary.maxwell.magnetic_field(RT_space, RT_space, NC_space, 1j*s)

	identity2=bempp.api.operators.boundary.sparse.identity(RT_space, RT_space, RT_space)

	identity= -bempp.api.operators.boundary.sparse.identity(RT_space, RT_space, NC_space)
	dof=NC_space.global_dof_count
	
	trace_fun= bempp.api.GridFunction(RT_space, coefficients=b[0:dof],dual_space=RT_space)

	zero_fun= bempp.api.GridFunction(RT_space,coefficients = b[dof:],dual_space=RT_space)
	
	#rhs=[trace_fun,zero_fun]
	id_discrete=identity2.weak_form()
	b[0:dof]=id_discrete*b[0:dof]
	blocks=np.array([[None,None], [None,None]])

	blocks[0,0] = -elec.weak_form()+0.1*s**0.5*identity2.weak_form()
	blocks[0,1] =  magn.weak_form()-1.0/2*identity.weak_form()
	blocks[1,0] = -magn.weak_form()-1.0/2*identity.weak_form()
	blocks[1,1] = -elec.weak_form()


	blocks_discrete=bempp.api.BlockedDiscreteOperator(blocks)
######## Saving the condition number and frequency : 
	
	from scipy.sparse.linalg import gmres
	lambda_data,info= gmres(blocks_discrete, b)
	#cond=np.linalg.cond(bempp.api.as_matrix(blocks_discrete))
	print("System solved !")
	import scipy.io,time
	mat_contents =scipy.io.loadmat('data/cond.mat')
	freqCond_old = mat_contents['freqCond']
	if s in freqCond_old[0]:
		print("Frequency already calculated")
	else:	
		tp0=time.time()
		blocks_mat=bempp.api.as_matrix(blocks_discrete)
	#	tp4=time.time()
		sigmas = scipy.linalg.svdvals(blocks_mat)
		norminv = min(sigmas)**(-1)
		normA = max(sigmas)
		cond=normA*norminv
		print("Freq: ",s ," Cond: ",cond)
	#	print(freqCond_old)
		freqCond=np.concatenate((freqCond_old,np.array([[s],[cond],[normA],[norminv]])),axis=1)
		scipy.io.savemat('data/cond.mat',dict(freqCond=freqCond))
#####################################################	
#print(np.linalg.norm(lambda_data))
	#print("I survived!")
	#from bempp.api.linalg import lu
	#lambda_data = lu(elec, trace_fun)
	#lambda_data.plot()
	#print("Norm lambda_data : ",np.linalg.norm(lambda_data))
	#if (np.linalg.norm(lambda_data)<10**-10):
	phigrid=bempp.api.GridFunction(RT_space,coefficients=lambda_data[0:dof],dual_space=RT_space)
	psigrid=bempp.api.GridFunction(RT_space,coefficients=lambda_data[dof:2*dof],dual_space=RT_space)


######## Create Points

#
#	x_a=-0.75
#        x_b=0.75
#        y_a=-0.25
#        y_b=1.25
##
#	x_a=-2
#        x_b=2
#        y_a=-2
#        y_b=2
#	n_grid_points=150
################################################
#        plot_grid = np.mgrid[y_a:y_b:1j*n_grid_points, x_a:x_b:1j*n_grid_points]
##       plot_grid = np.mgrid[-0.5:1:1j*n_grid_points, -1.5:1.5:1j*n_grid_points]
#        #print(plot_grid)
##        points = np.vstack( ( plot_grid[0].ravel() , plot_grid[1].ravel() , 0.25*np.ones(plot_grid[0].size) ) )
#
#	points = np.vstack( ( plot_grid[0].ravel()  , 0*np.ones(plot_grid[0].size) , plot_grid[1].ravel()) )
#point=np.array([[0],[0],[2]])
	slp_pot = bempp.api.operators.potential.maxwell.electric_field(RT_space, points, s*1j)

	dlp_pot = bempp.api.operators.potential.maxwell.magnetic_field(RT_space, points, s*1j)
	
	scattered_field_data = -slp_pot * phigrid+dlp_pot*psigrid
#	print("NORM COMBINED OPERATOR :" , np.linalg.norm(scattered_field_data)/np.linalg.norm(b))
#	print(scattered_field_data)
#	print("NORM ScatteredField :", np.linalg.norm(scattered_field_data))
#
#	print("s : ", s)
#	print("NORM B :" ,np.linalg.norm(b))
	if np.isnan(scattered_field_data).any():
		print("NAN Warning",s)
		print("NORM B :" ,np.linalg.norm(b))
		return np.zeros(n_grid_points**2*3)
	#print(scattered_field_data.reshape(3,1)[:,0])
	return scattered_field_data.reshape(n_grid_points**2*3,1)[:,0]

def scattering_solution(dx,N,T,m,points):
	import scipy.io
	import numpy as np
	mat_contents=scipy.io.loadmat('grids/TorusDOF896.mat')
	Nodes=np.array(mat_contents['Nodes']).T
	rawElements=mat_contents['Elements']
	for j in range(len(rawElements)):
		betw=rawElements[j][0]
		rawElements[j][0]=rawElements[j][1]
		rawElements[j][1]=betw
	Elements=np.array(rawElements).T
	Elements=Elements-1
	grid=bempp.api.grid_from_element_data(Nodes,Elements)
	
#	def tangential_trace(x, n, domain_index, result):
#		result[:] = n[1]
#
#	P1_space = bempp.api.function_space(grid,"P",1)
#	normal_fun = bempp.api.GridFunction(P1_space, fun=tangential_trace,dual_space=P1_space)
	#normal_fun.plot()

	#grid.plot()
	#grid=bempp.api.shapes.sphere(h=dx)
	rhs=create_rhs(grid,N,T,m)
	def ellipticSystem(s,b):
		return harmonic_calderon(s,b,grid,points)
	ScatOperator=Conv_Operator(ellipticSystem)
	#num_sol=ScatOperator.apply_convol(rhs,T)
	if (m==2):
		num_solStages=ScatOperator.apply_RKconvol(rhs,T,cutoff=10**(-16),method="RadauIIA-2")
	if (m==3):
		num_solStages=ScatOperator.apply_RKconvol(rhs,T,cutoff=10**(-16),method="RadauIIA-3")
	num_sol=np.zeros((len(num_solStages[:,0]),N+1))	
	num_sol[:,1:N+1]=np.real(num_solStages[:,m-1:N*m:m])
	return num_sol
import time

N=100
T=4

#Generate points

x_a=-1.5
x_b=1.5
y_a=-1.5
y_b=1.5
n_grid_points=300
nx=n_grid_points
nz=n_grid_points
#######################################
#Initialize empty file, which will be overwritten continously with condition numbers#and the frequencies
freqCond=1j*np.array([[0],[0],[0],[0]])
import scipy.io
scipy.io.savemat('data/cond.mat',dict(freqCond=freqCond))
###############
plot_grid = np.mgrid[y_a:y_b:1j*n_grid_points, x_a:x_b:1j*n_grid_points]
#plot_grid = np.mgrid[-0.5:1:1j*n_grid_points, -1.5:1.5:1j*n_grid_points]
#print(plot_grid)
points = np.vstack( ( plot_grid[0].ravel()  , 0*np.ones(plot_grid[0].size) , plot_grid[1].ravel()) )
	
evals=scattering_solution(1,N,T,3,points)
u_ges=np.zeros((n_grid_points**2,N+1))
for j in range(N+1):	
	#matplotlib inline
	import matplotlib
	from matplotlib import pylab as plt
	# Adjust the figure size in IPython
	matplotlib.rcParams['figure.figsize'] = (10.0, 8.0) 
	t=j*T*1.0/N

	def incident_field(x):
		return np.array([np.exp(-50*(x[2]-t+1)**2), 0. * x[2], 0. * x[2]])

	incident_field_data = incident_field(points)
	#scat_eval=np.zeros(nx*nz*3)
	#incident_field_data[radius<1]=np.nan
	scat_eval=evals[:,j].reshape(3,nx*nz)
#	print(scat_eval)
	field_data = scat_eval + incident_field_data
#	field_data = scat_eval 
#	field_data = incident_field_data 
#	print("Points: ")
#	print(points)
#	print("Data: ")
	#print(field_data)
	squared_field_density = np.sum(field_data * field_data,axis = 0)
	u_ges[:,j]=squared_field_density.T
	#squared_field_density=field_data[2,:]
	
	#squared_field_density[radius<1]=np.nan
	#print("MAX FIELD DATA: " , max(squared_field_density))	
	print(max(np.abs(squared_field_density)))
	
	plt.imshow(squared_field_density.reshape((nx, nz)).T,
	           cmap='coolwarm', origin='lower',
	           extent=[x_a, x_b, y_a, y_b])
	if j==10:
		plt.colorbar()
	plt.clim((-1,1))
	#plt.title("Squared Electric Field Density")

	#plt.savefig("data/wave_images/Screen_n{}.png".format(j))
import scipy.io
scipy.io.savemat('data/delta01_dof896.mat',dict(u_ges=u_ges,N=N,T=T,plot_grid=plot_grid,points=points))







#
#	tp1=time.time()
#	print("Dense Matrix assembled, time : ", tp1-tp0)
#
#	normA=np.linalg.norm(blocks_mat,ord=2)
#	tp2=time.time()
#	print("Norm A calculated, value ",normA,  " time : ",tp2-tp1)
#	cond=np.linalg.cond(blocks_mat)
#	tp3=time.time()	
#	print("Cond A calculated, value ", cond," time : ",tp3-tp2)
#	norminv=np.linalg.norm(np.linalg.inv(blocks_mat),ord=2)
#	print("Inv A calculated, direct : ", norminv, " Previous estimate : ", cond/normA)

