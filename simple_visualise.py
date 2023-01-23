import matplotlib.pyplot as plt
import numpy as np

pathname='sim_data/'
filename='init20230123_154621.dat'

with open(pathname + filename, "r") as f:
    data = np.loadtxt(f)
    x=data[:,1]
    y=data[:,2]
plt.figure()
plt.plot(x,y,'ro')
plt.show()