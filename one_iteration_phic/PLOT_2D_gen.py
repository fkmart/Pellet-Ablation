import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp, rc, pel_pot, lp, sig, t_start
import gauss_test_pot as gtp

direc = os.getcwd()
sub_dir = '/analysed_outputs/'
r_int = np.linspace(rp[t_start], rc[t_start], 602)

mid = 301

fig, ax  = plt.subplots()
ax2 = ax.twinx()
y = 1.00
for p in range(1, lp):
    pot = gtp.gauss_func(pel_pot[p],1.00, r_int[mid], r_int)
    file = np.loadtxt(direc + sub_dir + 'density_pot_test' + str(t_start) + 'pot' +str(int(pel_pot[p])) + 'sig_' + y + '.txt')
    ax.plot(file[:,0], file[:,1])
plt.show()