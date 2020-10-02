import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import *

mydir = './static_outputs_phi'

title = 'impacting_dist_fraction'

save_path = os.path.join(os.path.expanduser('~'), 'Pictures/t_static')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))

t = t_static
thing = np.loadtxt(os.path.join(mydir, 'retarded_flux_t'+str(t)+'.txt'))

fig, ax = plt.subplots() 
ax.plot(np.abs(thing[:,1]), thing[:,0])
ax.set_xlabel(r'$\phi_{\mathrm{pel}}/\mathrm{V}$', fontsize = 12)
ax.set_ylabel('Fraction of particles reaching pellet')
ax.set_yscale('log')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title(r'Fraction of distribution reaching pellet at $\tilde{t} = $' + str(t_static*dt) + ' for various retarding potentials')
plt.savefig(save_path+'/'+title+'.png', format = 'png', dpi = 1600)
plt.show()