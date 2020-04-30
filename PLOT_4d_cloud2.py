import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from gen_var import sig, pel_pot, rc, rp, t_start, r, lp 
import scipy.interpolate as spint
import os 
from matplotlib import cm

"""Challenges to getting a 4d plot made

1) Density/Bins on the colour axis 

2) Radius and Sigma on x and y

3) Potential on z axis but use flat layer like Alasdair's atmosphere plots
    - This negates the need to interpolate to a common grid for all combinations.
"""
t_start = 50
direc = os.getcwd()
load_dir = os.path.join(direc,'one_iteration_phic', 'analysed_outputs') + os.sep
save_dir = os.path.join(direc, 'Local_figures') + os.sep
fig = plt.figure(figsize = (10.0,7.5))
ax = fig.gca(projection='3d')
colors = ['dodgerblue','seagreen', 'gold', 'lilac']

sig = [0.5,0.6,0.7,0.8,0.9,1.0,1.25]
low = next(p[0] for p in enumerate(r) if p[1] > rp[t_start])
up = next(p[0] for p in enumerate(r) if p[1] > rc[t_start])
mid = int((up-low)/2)
r_internal = r[low:up]
bin_colour = np.zeros((len(sig),r_internal.size, lp))
X,Y = np.meshgrid(r_internal,sig)
Z = np.ones(bin_colour.size)
Z = np.reshape(Z,bin_colour.shape)
color_dimension = np.zeros((len(sig), r_internal.size,3))
for i in range(0, lp, 5):
    for j in range(0, len(sig)):
        Z[j,:,i] = pel_pot[i]/1000     
        file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start)+'pot' + str(pel_pot[i])+ 'sig'+str(sig[j]) + '.txt')
        f0,f1 = file[:,0], file[:,1]
        int_func = spint.interp1d(f0,f1, kind = 'cubic', fill_value = 'extrapolate')
        g = spint.PchipInterpolator(np.flip(f0, axis = 0), np.flip(f1, axis = 0), extrapolate = 'true')
        f1_full = int_func(r_internal)
        f1_full = g(r_internal)
        f1_full *= sum(file[:,1])/sum(f1_full) # normalised
        bin_colour[j,:,i] = f1_full
    color_dimension[:,:,int(i/5)] = bin_colour[:,:,i]
col = np.reshape(color_dimension,(X.size,3),'F')
vmi = np.min(color_dimension)
vma = np.max(color_dimension)

norm = matplotlib.colors.LogNorm(vmin = vmi,vmax = vma)

xs = X.shape
ys = Y.shape
zs = Z.shape

#Now produce the arrays for the "contours" that mark order of magnitude change
oom = 10**(np.linspace(-5,-3,3))
ch = np.zeros((len(sig),3,oom.size))
ch_val = np.zeros((len(sig),3,oom.size))
cent_line = np.zeros(len(sig))
cent_line[:] = r_internal[mid]
for k in range(0, oom.size):
    for i in range(0, 3):
        for j in range(0, len(sig)):
            ch[j,i,k] = next(p[0] for p in enumerate(color_dimension[j,:,i]) if p[1] > oom[k])    
            ch_val[j,i,k] = r_internal[int(ch[j,i,k])]
for i in range(0, 3):
    ax.plot_surface(X,Y,Z[:,:,-1 - 5*i], ccount = 30,
     facecolors = cm.magma(norm(color_dimension[:,:,-i])))
    if i==0:
        for b in range(0,oom.size):
            s = oom[b]
            ax.plot(ch_val[:,-i-1,b], sig[:], pel_pot[-1 - 5*i]*np.ones(len(sig))/1000, 
             color = colors[b], label = r'${:.0e}$'. format(s))
    else:
        for b in range(0,oom.size):
            ax.plot(ch_val[:,-i-1,b], sig[:], pel_pot[-1 - 5*i]*np.ones(len(sig))/1000, 
             color = colors[b])
    ax.plot(cent_line[:], sig, pel_pot[-1-5*i]*np.ones(len(sig))/1000, color = 'white', linestyle = '--',
     alpha = 0.6)

ax.set_xlabel(r'$\tilde{r}$', rotation = -0)
ax.set_ylabel(r'$\tilde{\sigma}$', rotation  = 90)
ax.set_zlabel(r'$V_{\mathrm{min}}/T_{MB}$', labelpad = 10, rotation = 0)
ax.set_xticklabels([1,2,3,4,5,6], rotation = 0, ha = 'left',va = 'bottom')
ax.set_yticklabels([0.50,0.75,1.00, 1.25], rotation = 0, ha = 'center', va = 'bottom')
ax.set_zticklabels([0,-1, -2], ha = 'center', va = 'center')
ax.set_zticks([0,-1,-2])
ax.set_yticks([0.5,0.75,1.0,1.25])
ax.tick_params(axis='x', pad=5)
ax.tick_params(axis = 'y', pad =4)
ax.tick_params(axis = 'z', pad = 5)
m = cm.ScalarMappable(cmap=plt.cm.magma, norm=norm)
m.set_array([])
cb = plt.colorbar(m, pad = 0.1)
cb.set_label('Electron Bins', rotation = 0, labelpad = 30, y = 0.5)
ax.view_init(20,-60)
ax.text2D(0.05,0.15, r'$\leftarrow$ to pellet', transform = ax.transAxes, rotation = -13)
ax.text2D(0.50, 0.05, r'$\rightarrow$ to plasma', transform= ax.transAxes, rotation = -13)
ax.text2D(0.05,0.90, r'$\tilde{\rho}_c(\tilde{r} = \tilde{r}_p) = 0.01$', transform = ax.transAxes)
ax.text2D(0.05, 0.85, r'$\tilde{\rho}_c(\tilde{r} = \tilde{r}_c) \approx 4.89 \times 10^{-4} $', transform = ax.transAxes)
ax.text2D(0.05,0.95, r'$T_{MB} = 1\mathrm{keV}$', transform = ax.transAxes)
plt.legend(title = 'EEDF Bin Markers', loc = 'best')
plt.savefig(save_dir + '4d_atmos_non-dim_markers.png', format = 'png', dpi = 1200)
plt.show()