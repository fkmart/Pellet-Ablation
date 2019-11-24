import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D 
import os

title = 'subplot_2d_cloud'

path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))

th = np.load('azi_mesh_2D.npy')
rad = np.load('rad_mesh_3D.npy')
dens = np.load('dens_mesh_3D.npy')

az = th[:,0]
lt = len(rad[2])

ind1 = 2
ind2 = 2
number = ind1*ind2
fig, axes = plt.subplots(nrows = ind1, ncols = ind2, subplot_kw ={'projection':'polar'})
ytick = [5,15,25,35]
ax1 = plt.subplot(221, projection = 'polar', xticks = [], yticks = ytick)
ax2 = plt.subplot(222, projection = 'polar')
ax3 = plt.subplot(223, projection = 'polar')
ax4 = plt.subplot(224, projection = 'polar')
c = 1
t = np.asarray([59,99,139,179])


im = ax1.pcolormesh(th, rad[:,:,t[0]],  dens[:,:,t[0]], cmap = 'viridis', norm = 
             colors.LogNorm(vmin = np.amin(dens[:,:,t[-1]]), vmax = np.amax(dens[:,:,t[0]])))
ax1.set_yticks(ytick)
ax1.set_xticks([])
ax1.grid()
ax1.text(45,40, r'$\tilde{t} = 0.3$')

ax2.pcolormesh(th, rad[:,:,t[1]],  dens[:,:,t[1]], cmap = 'viridis', norm = 
             colors.LogNorm(vmin = np.amin(dens[:,:,t[-1]]), vmax = np.amax(dens[:,:,t[0]])))
ax2.set_yticks(ytick)
ax2.set_xticks([])
ax2.grid()
ax2.text(45,40, r'$\tilde{t} = 0.5$')

ax3.pcolormesh(th, rad[:,:,t[2]],  dens[:,:,t[2]], cmap = 'viridis', norm = 
             colors.LogNorm(vmin = np.amin(dens[:,:,t[-1]]), vmax = np.amax(dens[:,:,t[0]])))
ax3.set_yticks(ytick)
ax3.set_xticks([])
ax3.grid()
ax3.text(45,40, r'$\tilde{t} = 0.7$')


ax4.pcolormesh(th, rad[:,:,t[3]],  dens[:,:,t[3]], cmap = 'viridis', norm = 
             colors.LogNorm(vmin = np.amin(dens[:,:,t[-1]]), vmax = np.amax(dens[:,:,t[0]])))
ax4.set_yticks(ytick)
ax4.set_xticks([])
ax4.grid()
ax4.text(45,40, r'$\tilde{t} = 0.9$')
"""
c = 0
for ax in axes.flat:
    ax.pcolormesh(th, rad[:,:,t[c]],  dens[:,:,t[c]], cmap = 'viridis', norm = 
             colors.LogNorm(vmin = np.amin(dens[:,:,t[-1]]), vmax = np.amax(dens[:,:,t[0]])))
    ax.set_yticks(ytick)
    ax.set_xticks([])
    ax.grid()
    c+=1
"""
cb_axes= fig.add_axes([0.9,0.1,0.03,0.8])
cb = fig.colorbar(im, cax = cb_axes ,ax = im)
#cb = plt.colorbar(im)
cb.set_label(r'$\tilde{\rho}$', rotation = 0, fontsize = '12', labelpad = -15)

"""
for tchoice in range(1, number + 1):
    r = rad[0,:,t[tchoice-1] ]
    #ax = fig.add_subplot(ind1,ind2,tchoice)
    #ax = plt.subplot(ind1,ind2, tchoice)
    #print(ax)
    #plt.subplot(projection = 'polar')
    
    #ax.plt.pcolormesh(th, rad[:,:,t[tchoice-1]],  dens[:,:,t[tchoice-1]], cmap = 'plasma', norm = 
    #           colors.LogNorm(vmin = np.amin(dens[:,:,t[tchoice-1]]), vmax = np.amax(dens[:,:,t[tchoice-1]])))

    #ax.plot(az, rad[:,:,t[tchoice-1]], color = 'white', ls = 'none')
    #ax.gridOn('True')
    #ax.xticks([])
    
    #ax.yticks((5,15,25,35))
   
    #plt.clim(1e-2, 1e-5)
   
    #cb = plt.colorbar()
    #cb.set_label(r'$\tilde{\rho}$', rotation = 0, fontsize = '12')
    #c +=1
    plt.show()
    plt.close()"""
#first one
#plt.savefig(path +'/' + title + '.png', format = 'png', dpi = 1600)
plt.show()