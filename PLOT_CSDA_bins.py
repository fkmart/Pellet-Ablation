import numpy as np 
import matplotlib.pyplot as plt 
import os 

fig, ax = plt.subplots()
direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration', 'neutral','1000eV','analysed_outputs','t_11772') + os.sep
file = np.loadtxt(load_dir + 'potential_t11772.txt')

pos = file[1,:]
bins = file[0,:]

print(sum(bins))
 
ax.plot(pos,bins)
#plt.yscale('log')
plt.show()