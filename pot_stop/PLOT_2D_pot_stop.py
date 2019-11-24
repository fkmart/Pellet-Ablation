import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp, rc, pel_pot, lp, sig, t_start
import gauss_test_pot as gtp

direc = os.getcwd()
sub_dir = '/pot_stop/analysed_data/'
r_int = np.linspace(rp[t_start], rc[t_start], 602)

mid = 301

fig, ax  = plt.subplots()

c = ['black', 'royalblue', 'chocolate', 'gold', 'darkorchid', 'limegreen','tomato', 'turquoise', 'midnightblue', 'firebrick', 'cadetblue']
y = format(1.00, '.2f')
for p in range(1, lp):
    pot = gtp.gauss_func(pel_pot[p],float(y), r_int[mid], r_int)
    x = str(int(pel_pot[p]))
    file = np.loadtxt(direc + sub_dir + 'elec_bins_peak' + x + 'sig_' + y + '.txt')
    ax.plot(file[:,0], file[:,1], color = c[p], label = r'$\phi_{\mathrm{max}} = $' + str(pel_pot[p]) + 'V')
plt.legend(ncol = 2, loc = 7)

#ax.set_yscale('log')
ax2 = ax.twinx()
for p in range(0, lp):
    pot = gtp.gauss_func(pel_pot[p],1.00, r_int[mid], r_int)
    ax2.plot(r_int, pot, color = c[p], linestyle = '--')

plt.show()