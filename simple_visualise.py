import matplotlib.pyplot as plt
import numpy as np

pathname='Sim_Data/'
filename='init.dat'

with open(pathname + filename, "r") as f:
    data = np.loadtxt(f)
    x=data[:,1]
    y=data[:,2]
plt.figure()
plt.plot(x,y,'ro')
plt.show()