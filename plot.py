import numpy as np
import matplotlib.pyplot as plt
from sys import exit


gr = np.loadtxt("data/gr.dat")

r = gr[:,0]
y = gr[:,1]

plt.plot(r,y)
plt.show()
