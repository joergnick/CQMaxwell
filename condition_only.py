import bempp.api
import numpy as np
from resource import *


import scipy.io
mat_contents=scipy.io.loadmat('data/cond_old.mat')
freqCond=mat_contents['freqCond']
freqs=freqCond[0,1:]
conds_old=freqCond[1,1:]
print("COND :", conds_old[4])
s=freqs[4]
import scipy.sparse
#
#def save_sparse_csr(filename,array):
#	np.savez(filename,data=array.data, indices=array.indices, indptr=array.indptr,shape=array.shape)
#def load_sparse_csr(filename):
#	loader = np.load(filename)	
#	return scipy.sparse.csr_matrix((loader['data'],loader['indices'],loader['indptr']),shape=loader['shape'])
#	
#Loading grid and reordering of Node count, to guarantee unified normal vector
#mat_contents=scipy.io.loadmat('grids/TorusDOF3392.mat')
mat_contents=scipy.io.loadmat('grids/TorusDOF896.mat')
#mat_contents=scipy.io.loadmat('grids/TorusDOF294.mat')
Nodes=np.array(mat_contents['Nodes']).T
rawElements=mat_contents['Elements']
for j in range(len(rawElements)):
	betw=rawElements[j][0]
	rawElements[j][0]=rawElements[j][1]
	rawElements[j][1]=betw
Elements=np.array(rawElements).T
Elements=Elements-1
grid=bempp.api.grid_from_element_data(Nodes,Elements)
	
##Bempp parameters

OrderQF = 7
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
#Bempp objects
NC_space=bempp.api.function_space(grid, "NC",0)
RT_space=bempp.api.function_space(grid, "RT",0)
dof=RT_space.global_dof_count

elec = -bempp.api.operators.boundary.maxwell.electric_field(RT_space, RT_space, NC_space,1j*s)

magn = -bempp.api.operators.boundary.maxwell.magnetic_field(RT_space, RT_space, NC_space, 1j*s)

identity2=bempp.api.operators.boundary.sparse.identity(RT_space, RT_space, RT_space)

identity= -bempp.api.operators.boundary.sparse.identity(RT_space, RT_space, NC_space)
dof=NC_space.global_dof_count

#blocks=np.array([[None,None], [None,None]])
#
#blocks[0,0] = -elec.weak_form()+0.1*s**0.5*identity2.weak_form()
#blocks[0,1] =  magn.weak_form()-1.0/2*identity.weak_form()
#blocks[1,0] = -magn.weak_form()-1.0/2*identity.weak_form()
#blocks[1,1] = -elec.weak_form()
#
#blocks_discrete=bempp.api.BlockedDiscreteOperator(blocks)
#
import time
#tp0=time.time()
#cond_old=np.linalg.cond(bempp.api.as_matrix(blocks_discrete))
#
tp1=time.time()
#print("Cond old : ",cond_old, " Time : ", tp1-tp0)
#
#save_sparse_csr('Mass.npz',id_mat)
#mass_mat=load_sparse_csr('Mass.npz')

elec_mat=bempp.api.as_matrix(elec.weak_form())
print(scipy.sparse.issparse(elec_mat))
import matplotlib.pyplot as plt
#plt.matshow(np.abs(elec_mat))
#plt.show()
tp2=time.time()
print("V form : ",elec_mat.shape, " Time needed : ", tp2-tp1)
magn_mat=bempp.api.as_matrix(magn.weak_form())

#scipy.sparse.save_npz('/data/magnetic.npz',magn_mat)
tp3=time.time()
print("K_form : ", magn_mat.shape," Time needed : ", tp3-tp2)
id_mat=bempp.api.as_matrix(identity.weak_form())
tp4=time.time()
print("Mass 1 : ", id_mat.shape," Time needed : ", tp4-tp3)
id_mat2=bempp.api.as_matrix(identity2.weak_form())
tp5=time.time()
print("Mass 2 : ", id_mat2.shape," Time needed : ", tp5-tp4)


#blockwise=np.array([[None,None], [None,None]])
#upleft=-elec_mat+0.1*s**0.5*id_mat2
#blockwise[0,0]= upleft
#upright=magn_mat-1.0/2*id_mat
#blockwise[0,1]= upright
#downleft=-magn_mat-1.0/2*id_mat
#blockwise[1,0]=downleft
#blockwise[1,1]=-elec_mat

#blockwise=np.array([[None,None], [None,None]])
blockwise=np.array([[-elec_mat+0.1*s**0.5*id_mat2,magn_mat-1.0/2*id_mat],[-magn_mat-1.0/2*id_mat,None]])

#blockwise[0,0]=-elec_mat+0.1*s**0.5*id_mat2

#blockwise[0,1]= magn_mat-1.0/2*id_mat

#blockwise[1,0]=-magn_mat-1.0/2*id_mat

blockwise[1,1]=-elec_mat
#bigmatrix=1j*np.zeros((2*dof,2*dof))
#
#bigmatrix[:dof,:dof] = -elec_mat+0.1*s**0.5*id_mat2
#bigmatrix[:dof,dof:] = magn_mat-1.0/2*id_mat
#bigmatrix[dof:2*dof,0:dof] = -magn_mat-1.0/2*id_mat
#bigmatrix[dof:2*dof,dof:2*dof]= -elec_mat

import scipy.sparse
#blocksparse = np.array([[None,None],[None,None]])
#blocksparse = np.array([[-elec_mat+0.1*s**0.5*id_mat2, magn_mat-1.0/2*id_mat],[-magn_mat-1.0/2*id_mat,-elec_mat]])

#completeBlockwise=scipy.sparse.bmat(blockwise).toarray()
completeBlockwise=scipy.sparse.csc_matrix(scipy.sparse.bmat(blockwise))

#A=scipy.sparse.csc_matrix([[1, 0, 0],[2,0,0],[1,1,1]],dtype=float)
#print("completeBlockwise : ", completeBlockwise)
sigmabig = scipy.sparse.linalg.svds(completeBlockwise, k=1,tol=0.0001,return_singular_vectors=False)
#conda=scipy.sparse.linalg.norm(completeBlockwise,ord=1)*scipy.sparse.linalg.norm(scipy.sparse.linalg.inv(completeBlockwise),ord=1)
#print(conda)
#x,istop,itn,normr,normar,norma,conda,normx = scipy.sparse.linalg.lsmr(completeBlockwise,2*1j*np.ones(2*dof))[:8]
#print(istop,itn,normr,norma,conda,np.linalg.norm(x),normx)
print("Success 1 : ", sigmabig)
#import inspect
#lines=inspect.getsource(scipy.sparse.linalg.svds)
#print(lines)
sigmasmall = scipy.sparse.linalg.svds(completeBlockwise, ncv=200,k=1,which='SM',tol=0.01,return_singular_vectors=False,maxiter=2000)
print(sigmabig*1.0/sigmasmall)
#print( " Success cond : ",sigmabig[-1],sigmabig[0], sigmabig[-1]*1.0/sigmabig[0])
#tpENDBLOCKWISE=time.time()
#print("Complete runtime for operator build blockwise: ", tpENDBLOCKWISE-tp1)
#
tp6=time.time()
print("Totaltime : ", tp6-tp1)

##print(dir(bempp.api.BlockedDiscreteOperator))
#import inspect
#lines=inspect.getsource(bempp.api.BlockedDiscreteOperator)
#print(lines)
#cond=np.linalg.cond(bempp.api.as_matrix(blocks_discrete))
#cond=np.linalg.cond(blocks_discrete)

#mat_contents =scipy.io.loadmat('data/cond.mat')
#freqCond_old = mat_contents['freqCond']
#
#print("Freq: ",s ," Cond: ",cond)
#print(freqCond_old)
#freqCond=np.concatenate((freqCond_old,np.array([[s],[cond]])),axis=1)
##	scipy.io.savemat('data/cond.mat',dict(freqCond=freqCond))
##
#
