import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = np.loadtxt("rhoz_dte-2.dat")
x = rho[:,0]
y = rho[:,1]

dx = x[1] - x[0]
N = dx * sum(y)
print(N)

plt.plot(x,y/N, label="ed bd e-2")


rho = np.loadtxt("rhoz.dat")
x = rho[:,0]
y = rho[:,1]

dx = x[1] - x[0]
N = dx * sum(y)
print(N)

plt.plot(x,y/N, label="ed bd e-3")

rho = np.loadtxt("../hard_sphere_mc/rhoz.dat")
x = rho[:,0]
y = rho[:,1] 

dx = x[1] - x[0]
N = dx * sum(y)
print(N)

plt.plot(x,y/N, label="mc")


plt.legend()
plt.show()
