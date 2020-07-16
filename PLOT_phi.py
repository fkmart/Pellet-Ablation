import numpy as np 
import matplotlib.pyplot as plt 
import os 

direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration', 'poisson','analysed_outputs','t_20') + os.sep
file = np.loadtxt(load_dir + 'potential_t20.txt')

fig, ax = plt.subplots() 
ax.plot(file[:,0], file[:,1])
plt.show()