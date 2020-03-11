import numpy as np
import os 
import matplotlib.pyplot as plt 

curdir = os.getcwd()
load_dir = os.path.join(curdir, 'one_iteration', 'analysed_outputs','t_40') + os.sep 
file = np.loadtxt(load_dir + 'density.txt')

fig, ax = plt.subplots() 
ax.plot(file[:,0], file[:,1])
#plt.savefig(load_dir + 't_50_dens_original.png', format = 'png', dpi = 1200)
plt.show()