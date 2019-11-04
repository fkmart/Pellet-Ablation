import numpy as np 
import os 
import matplotlib.pyplot as plt 
from gen_var import lp, t_static, pel_pot

mydir = './analysed_outputs'
title = 'min_imp_ener'

save_path = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))  

plot_file = []

for p in range(0, lp):
    file = np.loadtxt(mydir +  '/terminal_energy_t'+str(t_static)+'pot' + str(pel_pot[p]) +'.txt')
    plot_file.append(pel_pot[p])
    plot_file.append(file[0,0])

file = np.reshape(plot_file,(int(0.5*len(plot_file)), 2),'C')
fig, ax  = plt.subplots()
ax.plot(file[:,0], file[:,1])
#ax.set_yscale('log')
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.set_ylabel(r'$\tilde{\rho_e}$', fontsize = 12, rotation = 0)
ax.yaxis.set_label_coords(-0.07, 0.48)
ax.xaxis.set_label_coords(0.46, - 0.03)
#plt.title(r'Electron density in cloud at time $\tilde{t} = $' + str(t_static*dt) + ' for varying pellet potentials')
plt.legend()
#plt.savefig(save_path+'/'+title+'.png', format = 'png', dpi = 1600)
plt.show()