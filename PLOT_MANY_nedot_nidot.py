
import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp_hr, tf, delta_t, r0, N_0_sca, rpd_hr, many_start, tf

direc= os.getcwd() 

TYPE = ['neutral', 'poisson', 'sheath', 'ion']

limit = 2000
jump = 1
Nedot = np.zeros(limit)
Nidot = np.zeros(limit)
t_plot = np.linspace(many_start/10**6,(many_start + 50*2000)/10**6, limit)

fig, ax = plt.subplots()
plt.yscale('log')
plt.ylabel(r'$\dot{N}$', fontsize = 12, rotation = 0)
plt.xlabel(r'$\tilde{t}$')
ax.xaxis.set_label_coords(0.55,-0.03)
ax.yaxis.set_label_coords(-0.04,0.52)

load_dir = os.path.join(direc, 'many_iteration','ion','1000eV','analysed_outputs') + os.sep
for j in range(1, limit,jump):
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    N_e = float(file[3,4])
    N_i = float(file[3,5])
    Nedot[j] = N_e*N_0_sca/(delta_t*tf)    
    Nidot[j] = N_i*N_0_sca/(delta_t*tf)
ax.plot(t_plot, Nedot, label = r'$\dot{N}_e$')
ax.plot(t_plot,Nidot, label = r'$\dot{N}_i$')
plt.legend()
plt.savefig('electron_ion_ablation.png', format = 'png', dpi = 1400)
plt.show()