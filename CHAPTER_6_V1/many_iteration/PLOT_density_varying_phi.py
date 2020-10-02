import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import os
from gen_var import pel_pot, dt, t_static, lp 
import electron as el
mydir = './static_outputs_phi'

title = 'stopped_electron_density'

save_path = os.path.join(os.path.expanduser('~'), 'Pictures/t_static')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))

e_mid, e_bins, mb = el.dist_calc(el.e_dist, el.ener_res, el.e_bar)

fig, ax = plt.subplots()
for p in range(0,lp,10 ):
    thing2 = np.loadtxt(os.path.join(mydir,'terminal_energy_t'+str(t_static)+'pot'+str(pel_pot[p]) +'.txt'))
    thing3 = np.loadtxt(os.path.join(mydir,'stop_point_t' + str(t_static) +'pot'+str(pel_pot[p])+'.txt'))

    ind = np.where(thing3[:,1]==0)
    x = ind[0]
    x = x[0]
    sum1 = np.sum(e_bins[x:])

    thing = np.loadtxt(os.path.join(mydir,'density_t' +str(t_static) +'pot'+str(pel_pot[p])+'.txt'))
    sum2 = np.sum(thing[:,1])
    check = sum1 + sum2
    print('This should equal 1 : ' + str(check))
    ax.plot(thing[:,0], thing[:,1], label = r'$\phi_{\mathrm{pel}} = $' + str(pel_pot[p]))

#ax.set_yscale('log')
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.set_ylabel(r'$\tilde{\rho_e}$', fontsize = 12, rotation = 0)
ax.yaxis.set_label_coords(-0.07, 0.45)
ax.xaxis.set_label_coords(0.46, - 0.03)
plt.title(r'Electron density in cloud at time $\tilde{t} = $' + str(t_static*dt) + ' for varying pellet potentials')
plt.legend()
plt.savefig(save_path+'/'+title+'.png', format = 'png', dpi = 1600)
plt.show()

