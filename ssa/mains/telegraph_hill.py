from ssa._ssa import telegraph_constant
from  master import TwoStateConstMasterMatrixExp
from ode import TelegraphODE
import matplotlib.pyplot as plt
import time
import numpy as np

def get_discrete_values(times, x1, x2, x3, Nt, T=5):
    discrete_times = np.linspace(0, T, Nt)
    discrete_values = np.zeros((3, Nt))
    x3 = np.trim_zeros(x3,trim='b')
    nnonzero = len(x3)
    x1 = x1[:nnonzero]
    x2 = x2[:nnonzero]
    times = times[:nnonzero]
    if len(x1) > 0 and len(x2) > 0 and len(x3) > 0:
        discrete_values[0,0] = x1[0]
        discrete_values[1,0] = x2[0]
        discrete_values[2,0] = x3[0] 
        for i in range(1, Nt):
            dt = discrete_times[i] - discrete_times[i-1]
            time_idx = np.argwhere((times >= discrete_times[i-1]) & (times < discrete_times[i]))
            if np.any(time_idx):
                time_idx = np.squeeze(time_idx)
                discrete_values[0, i] = x1[time_idx]
                discrete_values[1, i] = x2[time_idx] 
                discrete_values[2, i] = x3[time_idx]
            else:
                discrete_values[0, i] = discrete_values[0, i-1]
                discrete_values[1, i] = discrete_values[1, i-1] 
                discrete_values[2, i] = discrete_values[2, i-1]
    else:
        discrete_values[0,:] = 1
        discrete_values[1,:] = 0
        discrete_values[2,:] = 0                     
    return discrete_times, discrete_values

def hill(x,K,n,beta):
    return beta/(1+(x/K)**n)

end_time = 48.0
k_on = 3.0
k_off = 3.0
ksyn = 5.0
kdeg = 3.0
K = 20
npow = 10
Nreps = 100
dt = 0.01
Nt = 2000
Nsamp = 50000

X = np.zeros((Nreps,3,Nsamp))
n = 0
while n < Nreps:
    print(f'Simulating Rep {n}')
    x1, x2, x3, times = telegraph_constant([end_time,k_on,k_off,ksyn,kdeg,Nt,K,npow])
    values = [x1, x2, x3]
    #print(times)
    try:
        discrete_times, discrete_values = get_discrete_values(times, x1, x2, x3, Nsamp,T=12.0)
        X[n,:,:] = discrete_values
        n += 1
    except: 
        print(f'Skipping rep {n} for reaction spacing')


pct_on = np.mean(X[:,1,:],axis=0)
pct_off = np.mean(X[:,0,:],axis=0)

fig, ax = plt.subplots(1,3,figsize=(10,3))
for n in range(Nreps):
    ax[0].plot(discrete_times,X[n,2,:])
ax[1].plot(discrete_times,np.mean(X[:,2,:],axis=0))
x = np.linspace(0,30,1000)
ax[2].plot(x,hill(x,K,npow,1))
ax[2].plot(x,hill(x,K,-npow,1))
plt.tight_layout()
plt.show()







