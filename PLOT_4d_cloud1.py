import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from gen_var import sig, pel_pot, rc, rp, t_start, r, lp 
import scipy.interpolate as spint
import os 

"""Challenges to getting a 4d plot made

1) Density/Bins on the colour axis 

2) Radius and Sigma on x and y

3) Potential on z axis but use flat layer like Alasdair's atmosphere plots
    - This negates the need to interpolate to a common grid for all combinations.
"""
direc = os.getcwd()
load_dir = os.path.join(direc,'one_iteration_phic', 'analysed_outputs') + os.sep
fig = plt.figure()
ax = fig.gca(projection='3d')

y = sig 
low = next(p[0] for p in enumerate(r) if p[1] > rp[t_start])
up = next(p[0] for p in enumerate(r) if p[1] > rc[t_start])
r_internal = r[low:up]
bin_colour = np.zeros((len(sig),r_internal.size))
X,Y = np.meshgrid(r_internal,sig)
for i in range(0, lp):
    for j in range(0, len(sig)):
        file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start)+'pot' + str(pel_pot[i])+ 'sig'+str(sig[j]) + '.txt')
        f0,f1 = file[:,0], file[:,1]
        int_func = spint.interp1d(f0,f1, kind = 'cubic', fill_value = 'extrapolate')
        f1_full = int_func(r_internal)
        f1_full /= sum(f1_full) # normalised
        bin_colour[j,:] = f1_full
    Z = pel_pot[i]*np.ones(X.size)
    Z = np.reshape(Z,X.shape)
    color_dimension = bin_colour # change to desired fourth dimension
    minn, maxx = color_dimension.min(), color_dimension.max()
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='plasma')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    ax.plot_surface(X,Y,Z,  facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
plt.show()