import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import expm

class TwoStateMaster:
    def __init__(self,r1,r2,a,b,k2t):
        self.r1 = r1
        self.r2 = r2
        self.a = a
        self.b = b
        self.k2t = k2t
    def W(self,t,r1,r2,a,b):
        k1t = a + b*(1-np.exp(-r1*t))*np.exp(-r2*t)
        k2t = self.k2t
        return np.array([[-k1t, k2t], [k1t, -k2t]])
    def dPdt(self,P,t,r1,r2,a,b):
        P = np.expand_dims(P,1)
        out = self.W(t,r1,r2,a,b) @ P
        return out
    def solve(self,P0,t):
        nt = len(t)
        dt = t[1]-t[0]
        P = np.zeros((2,nt))
        P[:,0] = P0
        for n in range(1,nt):
            eps = self.dPdt(P[:,n-1],t[n],self.r1,self.r2,self.a,self.b)*dt
            P[:,n] = P[:,n-1] + np.squeeze(eps)
        return P

class TwoStateConstMasterMatrixExp:
    def __init__(self,k_on,k_off):
        self.k_on = k_on
        self.k_off = k_off

    def solve(self,P0,t):
        k_on = self.k_on
        k_off = self.k_off
        W = np.array([[-k_on, k_off], [k_on, -k_off]])
        nt = len(t)
        P = np.zeros((2,nt))
        for n in range(nt):
            P[:,n] = P0 @ expm(W * t[n])
        return P

