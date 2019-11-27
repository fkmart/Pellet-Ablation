import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp, rc, pel_pot, lp, sig, t_start
import gauss_test_pot as gtp

direc = os.getcwd()
r_int = np.linspace(rp[t_start], rc[t_start], 602)

load_dir = os.path.join(direc, 'pot_stop', 'analysed_data') + os.sep
path = os.path.join(direc, 'pictures') + os.sep

mid = 301

fig, ax  = plt.subplots(figsize = (10.0,8.0))

c = ['black', 'royalblue', 'limegreen', 'gold', 'peru', 'tomato','midnightblue', 'forestgreen', 'y', 'saddlebrown', 'darkred']
y = format(1.00, '.2f')
for p in range(1, lp):
    pot = gtp.gauss_func(pel_pot[p],float(y), r_int[mid], r_int)
    x = str(int(pel_pot[p]))
    file = np.loadtxt(load_dir + 'elec_dens_peak' + x + 'sig_' + y + '.txt')
    ax.plot(file[:,1], file[:,0], color = c[p], label = str(pel_pot[p]) + 'V')
plt.legend(ncol = 1, loc = 'center right')

ax.set_yscale('log')
ax2 = ax.twinx()
for p in range(1, lp):
    pot = gtp.gauss_func(pel_pot[p],1.00, r_int[mid], r_int)
    ax2.plot(r_int, pot, color = c[p], linestyle = '--')

ax.set_xlabel(r'$r/r_0$', fontsize = 12)
ax.set_ylabel(r'$n_e / 10^{19}\mathrm{m}^{-3}$')
ax2.set_ylabel(r'$\phi$/V')
#plt.savefig(path + 'PES_dens_gauss_pot_many_sig.png', format = 'png', dpi = 800)
plt.show()