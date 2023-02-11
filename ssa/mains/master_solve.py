from ssa._ssa import telegraph
from  master import TwoStateMaster
import matplotlib.pyplot as plt
import time
import numpy as np

end_time = 24.0
k_off = 1.0
ksyn = 1
kdeg = 1
time_step = 0.005
nreps = 1000
dt = 0.1
x1_counts = []
x2_counts = []
x3_counts = []
reaction_times = []


a = 0
b = 4
r1 = 0.1
r2 = 0.1
Nt = 750
P0 = np.array([1,0])
t = np.linspace(0,end_time,int(end_time/dt))
tsm = TwoStateMaster(r1,r2,a,b,k_off)
P = tsm.solve(P0,t)


# Plot the average species counts versus time
fig, ax = plt.subplots(1,3,figsize=(9,2))
ax[0].plot(t,P[0,:],color='blue',linewidth=5.0,label=r'$P_{off}$')
ax[0].plot(t,P[1,:],color='red',linewidth=5.0,label=r'$P_{on}$')
ax[0].legend()
ax[0].set_xlabel('Time (hours)')
ax[0].set_ylabel('Probability')
plt.tight_layout()
plt.show()

