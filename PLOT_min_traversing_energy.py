import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import t_start, t_end, inc 

direc = os.getcwd() 

t = np.arange(float(t_start),float(t_end), 20.0)
t[:] = t[:]/t_end
ener = []
for i in range(t_start, t_end, 20):
    load_dir = os.path.join(direc, 'one_iteration', 'analysed_outputs','t_' + str(i)) + os.sep
    load = np.loadtxt(load_dir + os.sep + 'terminal_energy_neutral.txt')
    ener = np.append(ener, load[0,0])

fig, ax = plt.subplots() 
ax.plot(t,ener)
plt.xlabel(r'$\tilde{t}$', fontsize = 12)
plt.ylabel(r'$\vartheta$/eV', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.49, -0.04)
ax.yaxis.set_label_coords(-0.10, 0.40)
#plt.savefig('minimum_traversal_energy.png', format = 'png', dpi = 1400)
plt.show()