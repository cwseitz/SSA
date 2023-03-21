from SSA._SSA import photoswitch
from bin_ssa import bin_ssa
import matplotlib.pyplot as plt
import numpy as np

dt = 0.001 #us
T = 10e4*dt #100 ms

k12 = 400.0e-3 #ms^-1
k23 = 30.0e-3
k34 = 5.0e-3
k21 = 80.0e-3
k31 = 4.0e-3
k41 = 0.15e-3

trials = 1000

X = np.zeros((trials,4,round(T/dt)))
for n in range(trials):
    x1, x2, x3, x4, times = photoswitch([T,k12,k23,k34,k41,k31,k21])
    t_bins, x1_binned, x2_binned, x3_binned, x4_binned = bin_ssa(times, x1, x2, x3, x4, dt, T)
    X[n,0,:] = x1_binned
    X[n,1,:] = x2_binned
    X[n,2,:] = x3_binned
    X[n,3,:] = x4_binned
    
x1avg = np.mean(X[:,0,:],axis=0)
x2avg = np.mean(X[:,1,:],axis=0)
x3avg = np.mean(X[:,2,:],axis=0)
x4avg = np.mean(X[:,3,:],axis=0)


plt.plot(t_bins,x1avg,label=r'$x_{0}$',color='red')
plt.plot(t_bins,x2avg,label=r'$x_{1}$',color='blue')
plt.plot(t_bins,x3avg,label=r'$x_{2}$',color='purple')
plt.plot(t_bins,x4avg,label=r'$x_{3}$',color='black')
plt.xlabel('Time (ms)')
plt.legend()
plt.show()
