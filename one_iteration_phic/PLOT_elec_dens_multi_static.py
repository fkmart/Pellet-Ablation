import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import *
import os 

mydir = './analysed_outputs' 
savedir = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures')

fig, ax = plt.subplots() 
for p in range(0, len(pel_pot)):
    fil = np.loadtxt(os.path.join(mydir + '/real_density_t' + str(t_static) +'pot' +str(pel_pot[p]) + '.txt'))
    ax.plot(fil[:,1], fil[:,0], label = r'$\tilde{t} = $' + str(round((t_static/lt),2)))
ax.set_yscale('log')
ax.set_ylabel(r'$\tilde{n}_e$', fontsize = 12, rotation =0)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.xaxis.set_label_coords(0.53,-0.04)
ax.yaxis.set_label_coords(-0.06, 0.50)
plt.legend()
#plt.savefig(savedir+'/norm_e_dens_multi_static.png', format = 'png', dpi = 1200)
plt.show()