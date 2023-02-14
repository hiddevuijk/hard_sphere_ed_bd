import numpy as np
import matplotlib.pyplot as plt
from sys import exit


data = np.loadtxt("rhoz.dat")

x = data[:,0]
y = data[:,1]

plt.plot(x,y)
plt.show()
