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
