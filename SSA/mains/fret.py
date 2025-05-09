from SSA._SSA import fret_sim
import SSA
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

T=1000
k23=0.5; k34=0.5; k21=0.5
dtimes=[]; ftimes=[]; atimes=[]
for _ in range(10000):
    x1,x2,x3,x4,times=fret_sim([T,0,k23,k34,0,0,k21])
    X=np.array([x1,x2,x3,x4]).T
    for i,ti in enumerate(times):
        if X[i,0]==1: dtimes.append(ti)
        elif X[i,2]==1: ftimes.append(ti)
        elif X[i,3]==1: atimes.append(ti)

fig,ax=plt.subplots(3,1,figsize=(3,6),sharex=True)
bins=np.linspace(0,10,101)
t=0.5*(bins[:-1]+bins[1:])

Q=np.array([[-(k21+k23), k21,    k23,   0],
            [         0,   0,      0,    0],
            [         0,   0,   -k34, k34],
            [         0,   0,      0,    0]])

pdf_d=np.array([expm(Q*ti)[0,0]*k21 for ti in t])
pdf_f=np.array([expm(Q*ti)[0,0]*k23 for ti in t])
pdf_a=np.array([expm(Q*ti)[0,2]*k34 for ti in t])

lam=k21+k23
br_d=k21/lam; br_f=k23/lam

ax[0].hist(dtimes,bins=bins,density=True,alpha=0.5,color='k',label='sim')
ax[0].plot(t,pdf_d/br_d,label='theory')
ax[1].hist(ftimes,bins=bins,density=True,alpha=0.5,color='b',label='sim')
ax[1].plot(t,pdf_f/br_f,label='theory')
ax[2].hist(atimes,bins=bins,density=True,alpha=0.5,color='r',label='sim')
ax[2].plot(t,pdf_a/br_f,label='theory')

for a in ax:
    a.legend()

plt.tight_layout()
plt.show()
