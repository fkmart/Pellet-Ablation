import numpy as np
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d as m3d 
import matplotlib.cm as cm
from gen_var import sig, pel_pot, lp, rc, rp, t_start, x_res
import os
import scipy.interpolate as spint

direc = os.getcwd() 
load_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs') + os.sep 

fig, ax = plt.subplots(figsize = (10.0,8.00))

sig = 1.0
r = np.arange(rp[t_start], rc[t_start], x_res)
#pot = format(pel_pot[p], '.1f')
bin_arr = []
t_start = 50
for p in range(0,lp):
    pot = format(pel_pot[p], '.1f')
    k = format(sig , '.1f')
    file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start) +'pot'+str(pot) +'sig' + k +'.txt')
    X,Y = np.meshgrid(r, pel_pot)
    bins = file[:,1]
    f = spint.interp1d(file[:,0], bins, kind = 'cubic', fill_value = 'extrapolate')
    bins_plot = f(np.flip(r, axis = 0))
    bin_arr = np.append(bin_arr, bins_plot)
bin_arr = np.reshape(bin_arr, (lp, int(len(bin_arr)/lp)), order = 'C')
im = ax.imshow(bin_arr, interpolation = 'bilinear', extent = [r[-1], r[0],pel_pot[0], pel_pot[-1]],
                 cmap = cm.viridis, origin  = 'lower', aspect = 0.0005)
fig.colorbar(im, ax = ax, shrink = 0.42, aspect = 7)
plt.show()