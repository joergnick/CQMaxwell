
import numpy as np
def freq_der(s,b):
	return s**1*np.exp(-1*s)*b
import math
from RKconv_op import *

def create_timepoints(c,N,T):
	m=len(c)
	time_points=np.zeros((1,m*N))
	for j in range(m):
		time_points[0,j:m*N:m]=c[j]*1.0/N*np.ones((1,N))+np.linspace(0,1-1.0/N,N)
	return T*time_points


def create_rhs(N,T,m):
	if (m==2):
		c_RK=np.array([1.0/3,1])
	if (m==3):
		c_RK=np.array([2.0/5-math.sqrt(6)/10,2.0/5+math.sqrt(6)/10,1])

	rhs=np.zeros((1,N*m))
	time_points=create_timepoints(c_RK,N,T)
	for j in range(m*N):
		t=time_points[0,j]
		rhs[0,j]=np.sin(t)**6
	return rhs


def deriv_solution(N,T,m):
	rhs=create_rhs(N,T,m)
	def ellipticSystem(s,b):
		return harmonic_calderon(s,b,grid)
	ScatOperator=Conv_Operator(freq_der)
	#num_sol=ScatOperator.apply_convol(rhs,T)
	if (m==2):
		num_solStages=ScatOperator.apply_RKconvol(rhs,T,cutoff=10**(-8),show_progress=False,method="RadauIIA-2")
	if (m==3):
		num_solStages=ScatOperator.apply_RKconvol(rhs,T,cutoff=10**(-8),show_progress=False,method="RadauIIA-3")
	num_sol=np.zeros((1,N+1))	
	num_sol[:,1:N+1]=np.real(num_solStages[m-1:N*m:m])
	return num_sol

import time


T=2
N_ref=2**12
tt_ref=np.linspace(0,T,N_ref+1)
m=3


sol_ref=deriv_solution(N_ref,T,m)
Am_time=8
#Am_space=1
#Am_time=8
tau_s=np.zeros(Am_time)
errors=np.zeros(Am_time)

m=2
for ixTime in range(Am_time):
	N=8*2**(ixTime)
	tau_s[ixTime]=T*1.0/N
	tt=np.linspace(0,T,N+1)

## Rescaling reference solution:		
	speed=N_ref/N
	resc_ref=np.zeros((3,N+1))
#	resc_ref=sol_ref
	for j in range(N+1):
		resc_ref[:,j]      = sol_ref[:,j*speed]
	#num_sol = calc_ref_sol(N,dx,F_transfer)	
	num_sol  = deriv_solution(N,T,m)
#	plt.plot(tt,num_sol[0,:]**2+num_sol[1,:]**2+num_sol[2,:]**2)
#	plt.plot(tt_ref,sol_ref[0,:]**2+sol_ref[1,:]**2+sol_ref[2,:]**2,linestyle='dashed')
	#plt.show()
	errors[ixTime]=np.max(np.abs(resc_ref-num_sol))
	print(errors)
	import scipy.io
#	scipy.io.savemat('data/Err_data_delta01.mat', dict( ERR=errors,h_s=h_s,tau_s=tau_s))
	#scipy.io.savemat('data/Err_data_delta0p1_long.mat', dict( ERR=errors,h_s=h_s,tau_s=tau_s))
#end=time.time()
#print("Script Runtime: "+str((end-start)/60) +" Min")
import matplotlib.pyplot as plt
plt.loglog(tau_s,errors)
plt.loglog(tau_s,tau_s**3,linestyle='dashed')
plt.loglog(tau_s,tau_s**2,linestyle='dashed')
plt.loglog(tau_s,tau_s**1,linestyle='dashed')
plt.show()
