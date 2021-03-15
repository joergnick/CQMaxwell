import scipy.io
import numpy as np

#v=1j*np.zeros((2,1))
#scipy.io.savemat('data/test.mat',dict(v=v))
mat_contents=scipy.io.loadmat('data/cond.mat')
freqCond=mat_contents['freqCond']
freqs = np.concatenate((freqCond[0,1:],np.conj(freqCond[0,1:])))
conds = np.concatenate((freqCond[1,1:],freqCond[1,1:]))
norms = np.concatenate((freqCond[2,1:],freqCond[2,1:]))
norminvs = np.concatenate((freqCond[3,1:],freqCond[3,1:]))

n=len(freqs)
posImagpart=[]
negImagpart=[]

for j in range(n):
	if np.imag(freqs[j])>=0:
		posImagpart.append(j)
	else:
		negImagpart.append(j)

freqsPos=freqs[posImagpart]
condsPos=conds[posImagpart]
normsPos=norms[posImagpart]
norminvsPos=normsinvs[posImagpart]

indicesPos=np.argsort(np.real(freqsPos))

freqsPosOrdered=freqsPos[indicesPos]
condsPosOrdered=condsPos[indicesPos]
normsPosOrdered=normsPos[indicesPos]
norminvsPosOrdered=norminvsPos[indicesPos]

import matplotlib.pyplot as plt
#plt.scatter(np.real(freqsPos),np.imag(freqsPos))
#plt.show()

plt.semilogy(condsPosOrdered,linestyle='-')
#plt.semilogy(np.real(freqsPos))
plt.show()
#v_old=mat_contents['v']
#n=len(v_old[0,:])
#x=1j
#y=1+1j
#v_old=[[],[]]
#v_new=np.concatenate((v_old,np.array([[x],[y]])),axis=1)
##v_new=np.zeros((2,n+1))
##v_new[:,:n]=v_old
##v_new[:,n]=np.array([x,y])

