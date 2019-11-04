import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import t_static, dt
import electron  


mydir = './analysed_outputs'

title = 'impacting_energy_flux'

save_path = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))
t = t_static 

thing = np.loadtxt(os.path.join(mydir, 'retarded_ener_flux'+str(t)+'.txt')) 

fig, ax = plt.subplots()
ax.plot(np.abs(thing[:,1]), thing[:,0])
ax.set_xlabel(r'$\phi_{\mathrm{pel}}$/V', fontsize = 12)
ax.set_ylabel(r'Energy flux arriving at pellet / $\mathrm{eVm}^{-2}\mathrm{s}^{-1}$')
ax.set_yscale('log')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.title(r'Energy flux at pellet surface at $\tilde{t} = $' + str(t_static*dt) + ' for various retarding potentials')
plt.savefig(save_path+'/'+title+'.png', format = 'png', dpi = 1600)
plt.show()