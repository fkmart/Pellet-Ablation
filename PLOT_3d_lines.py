import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
from gen_var import p_inc, pel_pot, lp, t_start, rp, rc, sig, r, n_r
import os 
import gauss_test_pot
direc = os.getcwd() 

load_dir = direc + '/one_iteration_phic/analysed_outputs/'
path = '/home/kyle/Pictures/'
i = t_start 

low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
up = next(p[0] for p in enumerate(r) if p[1] > rc[i])

r_internal = r[low:up]
mid = int((up-low)/2.0)

lines = []
fig = plt.figure()
ax = plt.axes(projection = '3d')
#ax2 = plt.axes(projection = '3d', sharex = 'true')
c = ['black', 'royalblue', 'chocolate', 'gold', 'darkorchid', 'limegreen','tomato', 'turquoise', 'midnightblue', 'firebrick', 'cadetblue']
for a in range(0,len(sig)):
    pot = gauss_test_pot.gauss_func(0.5, sig[a], r_internal[mid], r_internal)
    for p in range(0, lp, p_inc):
        file = np.loadtxt(load_dir + 'real_density_pot_test_t'+str(i) +'pot'+str(pel_pot[p])+'sig' + str(sig[a])+'.txt')
        x = file[:,0]
        y = file[:,1]
        z = np.ones(len(x)) * sig[a]
        z_pot = np.ones(len(r_internal)) * sig[a]
        if (a==0):
            lines += ax.plot3D(y,z,x, color = c[p],label = r'$\phi = $' + str(pel_pot[p]) +'V', linewidth = 1.0)
            ax.plot(r_internal,z_pot,pot, linestyle = '--', color = 'black')
        else:
            ax.plot3D(y,z,x, color = c[p],label = r'$\phi = $' + str(pel_pot[p]), linewidth = 1.0) 
            ax.plot(r_internal,z_pot,pot, linestyle = '--', color = 'black')
    labs = [l.get_label() for l in lines]
ax.set_xlabel(r'$r/r_0$')
ax.set_ylabel(r'$\sigma /r_0$', rotation = 0)
ax.set_zlabel(r'$n_e/10^{19}\mathrm{m}^{-3}$', rotation = 0)
#ax.set_zlabel('Electron bins')
ax.legend(lines, labs, loc = 9, ncol = 4)
ax.view_init(30,120)
plt.savefig(path + 'gauss_pot_elec_dens_wpot.png', format = 'png', dpi = 1400)
plt.show()        

"""
a = 3
pot = gauss_test_pot.gauss_func(1000.0, sig[a], r_internal[mid], r_internal)
fig, ax  = plt.subplots() 
for p in range(0, lp, p_inc):
    file = np.loadtxt(load_dir + 'density_pot_test_t'+str(i) +'pot'+str(pel_pot[p])+'sig' + str(sig[a])+'.txt')
    x = file[:,0]
    y = file[:,1]
    ax.plot(x,y, label = r'$\phi = $' + str(pel_pot[p]))
plt.legend()
ax2 = ax.twinx()
ax2.plot(r_internal, pot, '--')
plt.show()
"""