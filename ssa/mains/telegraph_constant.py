from ssa._ssa import telegraph_constant
from  master import TwoStateConstMasterMatrixExp
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


end_time = 5.0
k_on = 1.0
k_off = 1.0
ksyn = 3.0
kdeg = 1.0
time_step = 0.005
Nreps = 200
dt = 0.1
Nt = 500
Nsamp = 20000

P0 = np.array([1,0])
t = np.linspace(0,end_time,Nt)
tsm = TwoStateConstMasterMatrixExp(k_on,k_off)
P = tsm.solve(P0,t)

X = np.zeros((Nreps,3,Nsamp))
for n in range(Nreps):
    print(f'Simulating Rep {n}')
    x1, x2, x3, times = telegraph_constant([end_time,k_on,k_off,ksyn,kdeg,Nt])
    values = [x1, x2, x3]
    discrete_times, discrete_values = get_discrete_values(times, x1, x2, x3, Nsamp)
    X[n,:,:] = discrete_values


pct_on = np.mean(X[:,1,:],axis=0)
pct_off = np.mean(X[:,0,:],axis=0)

# Plot the average species counts versus time
fig, ax = plt.subplots(figsize=(6,2))
ax.plot(t,P[0,:],color='blue',linewidth=5.0,label=r'$P_{off}$')
ax.plot(t,P[1,:],color='red',linewidth=5.0,label=r'$P_{on}$')
ax.legend()
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Probability')
ax.scatter(discrete_times,pct_on,color='red')
ax.scatter(discrete_times,pct_off,color='blue')
plt.tight_layout()
plt.show()








