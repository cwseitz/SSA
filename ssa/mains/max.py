import numpy as np
import matplotlib.pyplot as plt

a = np.linspace(0, 5, 100)
b = np.linspace(0, 5, 100)
A, B = np.meshgrid(a, b)

t = -np.log(A / (A + B)) / (A + B)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(A, B, t, cmap='coolwarm')
plt.show()

