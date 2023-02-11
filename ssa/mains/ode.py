import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class TelegraphConstODE:

    def __init__(self,k1,k2,k3,k4):

        # Define rate constants
        self.k1 = k1 #switch on
        self.k2 = k2 #switch off
        self.k3 = k3 #transcription
        self.k4 = k4 #degradation


    def ode_system(self, t, X, k1, k2, k3, k4):
        x, y = X
        dx = k1 * (1 - x) - k2 * x
        dy = k3 * x - k4 * y
        return [dx, dy]
        
    def solve(self,t_span=(0,25),X0=[0.1, 0.1]):

        k1, k2, k3, k4 = self.k1, self.k2, self.k3, self.k4
        solution =\
        solve_ivp(lambda t, X: self.ode_system(t, X, k1, k2, k3, k4), t_span, X0)
        
        t = solution.t
        x = solution.y[0]
        y = solution.y[1]
        return t,x,y


class TelegraphHillODE:

    def __init__(self,k1,k2,k3,k4):

        # Define rate constants
        self.k1 = 0.5 #switch on
        self.k2 = 0.5 #switch off
        self.k3 = 10.0 #transcription
        self.k4 = 1.0 #degradation
        self.K = 5
        self.n = 3

    def hill(self,x,K,n):
        return 1/(1+(x/K)**n)

    # Define ODE system
    def ode_system(self, t, X, k1, k2, k3, k4):
        x, y = X
        dx = k1 * (1 - x) - k2 * x
        dy = k3 * x - k4 * y
        return [dx, dy]
        
    def solve(self):

        X0 = [0.1, 0.1]
        k1, k2, k3, k4 = self.k1, self.k2, self.k3, self.k4
        t_span = (0, 25)
        solution =\
        solve_ivp(lambda t, X: self.ode_system(t, X, k1, k2, k3, k4), t_span, X0)
        
        t = solution.t
        x = solution.y[0]
        y = solution.y[1]
        num = np.linspace(0,20,100)

        fig, ax = plt.subplots(1,2,figsize=(8,3))
        ax[0].plot(t,x)
        ax[1].plot(t,y)
        plt.tight_layout()
        plt.show()
