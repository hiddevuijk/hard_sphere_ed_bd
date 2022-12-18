import numpy as np
import matplotlib.pyplot as plt
from sys import exit


pcorr = np.loadtxt("data/pcorr.dat")
r = pcorr[:,0]
gr = pcorr[:,1]

plt.plot(r,gr)
plt.show()
