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
X,Y = np.meshgrid(r_internal,sig)

bin_colour = np.zeros(X.shape)
for i in range(10, lp, lp):
    Z = np.ones(X.size)
    Z = np.reshape(Z,X.shape)
    Z[:,:] *= pel_pot[i] 
    for j in range(0, len(sig)):   
        file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start)+'pot' + str(pel_pot[i])+ 'sig'+str(sig[j]) + '.txt')
        f0,f1 = file[:,0], file[:,1]
        int_func = spint.interp1d(f0,f1, kind = 'cubic', fill_value = 'extrapolate')
        g = spint.PchipInterpolator(np.flip(f0, axis = 0), np.flip(f1, axis = 0), extrapolate = 'true')
        f1_full = int_func(r_internal)
        f1_full = g(r_internal)
        #f1_full *= sum(file[:,1])/sum(f1_full) # normalised
        check = np.flip(f1_full, axis = 0)
        bin_colour[j,:] = f1_full
    color_dimension = np.copy(bin_colour) # change to desired fourth dimension
    #color_dimension = np.reshape(color_dimension,(len(sig)* r_internal.size,3), 'F')
    minn, maxx = color_dimension.min(), color_dimension.max()
    #norm = matplotlib.colors.Normalize(minn, maxx)
    norm = matplotlib.colors.LogNorm(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='viridis')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    ax.plot_surface(X,Y,Z,  facecolors=fcolors, vmin=minn, vmax=maxx, shade=False, ccount = r_internal.size)
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\tilde{\sigma}$')
ax.set_zlabel('Peak Potential/V')
fig.colorbar(m)
plt.show()