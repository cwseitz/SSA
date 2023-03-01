from ssa._ssa import photoswitch
from bin_ssa import bin_ssa
import matplotlib.pyplot as plt
import numpy as np

T = 100.0
k1 = 0.2
k2 = 0.1
k3 = 0.05
k4 = 0.05
dt = 0.01
trials = 1000
X = np.zeros((trials,3,round(T/dt)))
for n in range(trials):
    x1, x2, x3, times = photoswitch([T,k1,k2,k3,k4])
    t_bins, x1_binned, x2_binned, x3_binned = bin_ssa(times, x1, x2, x3, dt, T)
    X[n,0,:] = x1_binned
    X[n,1,:] = x2_binned
    X[n,2,:] = x3_binned

x1avg = np.mean(X[:,0,:],axis=0)
x2avg = np.mean(X[:,1,:],axis=0)
x3avg = np.mean(X[:,2,:],axis=0)

plt.plot(t_bins,x1avg)
plt.plot(t_bins,x2avg)
plt.plot(t_bins,x3avg)
plt.show()
