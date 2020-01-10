import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.interpolate as spint 
from matplotlib.colors import LogNorm
from gen_var import sig, pel_pot, rc, rp, t_start, r, lp 

direc = os.getcwd()
load_dir = os.path.join(direc,'one_iteration_phic', 'analysed_outputs') + os.sep 

y = sig 
low = next(p[0] for p in enumerate(r) if p[1] > rp[t_start])
up = next(p[0] for p in enumerate(r) if p[1] > rc[t_start])
r_internal = r[low:up]
X,Y = np.meshgrid(r_internal,sig)

bin_colour = np.zeros(X.shape) 

ind = -1

Z = np.ones(X.size)
Z = np.reshape(Z,X.shape)
Z[:,:] *= pel_pot[ind] 
fig, ax = plt.subplots()
for j in range(0, len(sig)): 
    file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start)+'pot' + str(pel_pot[ind])+ 'sig'+str(sig[j]) + '.txt')
    f0,f1 = file[:,0], file[:,1]
    #int_func = spint.interp1d(f0,f1, kind = 'cubic', fill_value = 'extrapolate')
    g = spint.PchipInterpolator(np.flip(f0, axis = 0), np.flip(f1, axis = 0), extrapolate = 'true')
    #f1_full = int_func(r_internal)
    f1_full = g(r_internal)
    #f1_full *= sum(file[:,1])/sum(f1_full) # normalised
    check = np.flip(f1_full, axis = 0)
    bin_colour[j,:] = f1_full
    
cax = ax.imshow(bin_colour, cmap = 'viridis', aspect = 100, norm = LogNorm(vmin = bin_colour.min(), vmax= bin_colour.max()), interpolation = 'none')
cbar = fig.colorbar(cax)
plt.show()