from ssa._ssa import two_state
from  master import TwoStateMaster
import matplotlib.pyplot as plt
import time
import numpy as np

end_time = 24.0
k_on = 1.0
k_off = 1.0
time_step = 0.005
nreps = 2000
dt = 0.1
x1_counts = []
x2_counts = []
reaction_times = []


a = 0
b = 4
r1 = 0.2
r2 = 0.2
P0 = np.array([1,0])
t = np.linspace(0,end_time,int(end_time/dt))
tsm = TwoStateMaster(r1,r2,a,b,k_off)
P = tsm.solve(P0,t)

for n in range(nreps):
    systime = time.time()
    x1, x2, times = two_state([end_time,time_step,k_on,k_off])
    x1_counts.append(x1)
    x2_counts.append(x2)
    reaction_times.append(times)
    
x1_counts = np.array(x1_counts)
x2_counts = np.array(x2_counts)
reaction_times = np.array(reaction_times)

# Calculate the number of time steps
num_time_steps = int(np.ceil(reaction_times.max() / dt))
time = np.arange(0, num_time_steps) * dt

# Initialize an array to store the average species counts
avg_x1_counts = np.zeros((num_time_steps,))
avg_x2_counts = np.zeros((num_time_steps,))

# Initialize an array to store the number of species counts at each time step
num_x1_counts = np.zeros((num_time_steps,))
num_x2_counts = np.zeros((num_time_steps,))

# Determine the bin indices for each reaction time
bin_indices = np.digitize(reaction_times.flatten(), np.arange(0, num_time_steps) * dt)

# Loop over all replicates and reactions
for i, j in zip(x1_counts.flatten(), bin_indices):
    # Add the species count to the running sum for this time step
    avg_x1_counts[j - 1] += i
    # Increment the number of species counts at this time step
    num_x1_counts[j - 1] += 1

# Loop over all replicates and reactions
for i, j in zip(x2_counts.flatten(), bin_indices):
    # Add the species count to the running sum for this time step
    avg_x2_counts[j - 1] += i
    # Increment the number of species counts at this time step
    num_x2_counts[j - 1] += 1

# Divide the running sum by the number of species counts to get the average
avg_x1_counts /= num_x1_counts
avg_x2_counts /= num_x2_counts

# Plot the average species counts versus time
import matplotlib.pyplot as plt
plt.plot(t,P[0,:],color='blue',linewidth=5.0,label=r'$P_{on}$')
plt.plot(t,P[1,:],color='red',linewidth=5.0,label=r'$P_{off}$')
plt.scatter(time[1:], avg_x1_counts[1:],s=7,color='blue')
plt.scatter(time[1:], avg_x2_counts[1:],s=7,color='red')
plt.legend()
plt.xlabel('Time (hours)')
plt.ylabel('Probability')
plt.show()






