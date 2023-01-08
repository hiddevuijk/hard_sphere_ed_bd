import numpy as np
import matplotlib.pyplot as plt
from sys import exit


rho = np.loadtxt("rhoz.dat")
x = rho[:,0]
y = rho[:,1]

dx = x[1] - x[0]
N = dx * sum(y)
print(N)

plt.plot(x,y, label="ed bd")

rho = np.loadtxt("../hard_sphere_mc/rhoz.dat")
x = rho[:,0]
y = rho[:,1] 

dx = x[1] - x[0]
N = dx * sum(y)
print(N)

plt.plot(x,y, label="mc")


plt.legend()
plt.show()
