import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import scipy.integrate as spinteg
import os
from gen_var import pel_pot, dt, t_static, lp 
import electron as el
mydir = './static_outputs_phi'

title = 'real_stopped_electron_density'

save_path = os.path.join(os.path.expanduser('~'), 'Pictures/t_static')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))

e_mid, e_bins, mb = el.dist_calc(el.e_dist, el.ener_res, el.e_bar)

fig, ax = plt.subplots()
for p in range(0,lp,10 ):
    thing = np.loadtxt(os.path.join(mydir,'real_density_t' +str(t_static) +'pot'+str(pel_pot[p])+'.txt'))
    check = spinteg.simps(thing[:,0], thing[:,1])
    print('This should equal 1 : ' + str(check))
    ax.plot(thing[:,1], thing[:,0], label = r'$\phi_{\mathrm{pel}} = $' + str(pel_pot[p]))

#ax.set_yscale('log')
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.set_ylabel(r'$\tilde{\rho_e}$', fontsize = 12, rotation = 0)
ax.yaxis.set_label_coords(-0.07, 0.48)
ax.xaxis.set_label_coords(0.46, - 0.03)
plt.title(r'Electron density in cloud at time $\tilde{t} = $' + str(t_static*dt) + ' for varying pellet potentials')
plt.legend()
plt.savefig(save_path+'/'+title+'.png', format = 'png', dpi = 1600)
plt.show()

