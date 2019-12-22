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
    file = np.loadtxt(load_dir + 'stop_point_peak' + x + 'sig_' + y + '.txt')
    ax.scatter(file[:,1], file[:,0], marker = 'x', color = c[p], label = str(pel_pot[p]) + 'V')
plt.legend(ncol = 1, loc = 'upper right')

ax.set_yscale('log')

ax.axvline(r_int[mid], linestyle = '--', color = 'black',label = 'Gaussian Centre')
ax.axvline(r_int[0], linestyle = ':', color = 'k')
ax.axvline(r_int[-1], linestyle = ':', color = 'k')
plt.text(r_int[mid + int(0.05*mid) ],10**4, 'Gaussian Centre')
plt.text(r_int[-20], 0.6*10**2, r'$r_c$', fontsize = 12)
plt.text(r_int[10], 10**2, r'$r_p$', fontsize = 12)

ax.set_xlabel(r'$r/r_0$', fontsize = 12)
ax.set_ylabel('Initial Electron Energy/eV')

plt.savefig(path + 'PES_stop_point_gauss_pot_many_sig.png', format = 'png', dpi = 800)
plt.show()