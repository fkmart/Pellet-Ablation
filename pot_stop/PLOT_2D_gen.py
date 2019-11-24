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
ax2 = ax.twinx()
for a in range(0, len(sig)):
    y = format(sig[a],'.2f')
    for p in range(1, lp):
        pot = gtp.gauss_func(pel_pot[p],sig[a], r_int[mid], r_int)
        file = np.loadtxt(direc + sub_dir + 'elec_bins_peak' +str(int(pel_pot[p])) + 'sig_' + y + '.txt')
        ax.plot(file[:,0], file[:,1])
plt.show()