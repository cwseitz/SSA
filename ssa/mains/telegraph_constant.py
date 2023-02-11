from ssa._ssa import telegraph_constant
from ode import TelegraphConstODE
from discrete import get_discrete_values
import matplotlib.pyplot as plt
import time
import numpy as np

##########################
# Params
##########################

end_time = 48.0
k1 = 2.0 #on
k2 = 2.0 #off
k3 = 5.0 #syn
k4 = 0.0 #deg
Nreps = 100
dt = 0.01
Nt = 2000
Nsamp = 50000

##########################
# Solve ODE system
##########################

ode = TelegraphConstODE(k1,k2,k3,k4)
odet,odex,odey = ode.solve(X0=[1,0])

##########################
# Simulate SSA
##########################


X = np.zeros((Nreps,3,Nsamp))
n = 0
while n < Nreps:
    print(f'Simulating Rep {n}')
    x1, x2, x3, times = telegraph_constant([end_time,k1,k2,k3,k4,Nt])
    values = [x1, x2, x3]
    try:
        discrete_times, discrete_values = get_discrete_values(times, x1, x2, x3, Nsamp,T=12.0)
        X[n,:,:] = discrete_values
        n += 1
    except: 
        print(f'Skipping rep {n} for reaction spacing')


pct_on = np.mean(X[:,1,:],axis=0)
pct_off = np.mean(X[:,0,:],axis=0)

fig, ax = plt.subplots(1,2,figsize=(6,3))
for n in range(Nreps):
    ax[0].plot(discrete_times,X[n,2,:])
ax[1].plot(discrete_times,np.mean(X[:,2,:],axis=0))
ax[1].plot(odet,odey,color='red')
x = np.linspace(0,30,1000)
plt.tight_layout()
plt.show()







