import numpy as np 
import matplotlib.pyplot as plt
import os 

no = '161'

c = os.getcwd() + os.sep
cd = os.path.join(c,'one_iteration','analysed_outputs','t_40')
file = np.loadtxt(cd + os.sep+'neutral_rk\density.txt')
file2 = np.loadtxt(cd +  os.sep+'sheath_rk\density.txt')
file3 = np.loadtxt(cd + os.sep + 'neutral_rkf\density.txt')
file4 = np.loadtxt(cd + os.sep + 'sheath_rkf\density.txt')

fig, ax = plt.subplots()
ax.plot(file[:,1], file[:,0])
ax.plot(file2[:,1], file2[:,0])
ax.plot(file3[:,1], file3[:,0])
ax.plot(file4[:,1], file4[:,0])
plt.yscale('log')
plt.show()