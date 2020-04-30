import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp, rc 
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D 

fig = plt.figure() 
ax = fig.add_subplot(111, polar = 'true')

i = 300
rc_par = rc[i] # analytically known
rp_par = rp[i] # analytically known

r = np.linspace(0.0, np.round(rc[i]), num = 512, endpoint = 'true')
ind_cpar = (np.abs(r - rc_par)).argmin()
ind_ppar = (np.abs(r - rp_par)).argmin()
theta = np.linspace(0.0, 2.0*np.pi, num = 512)

R,T = np.meshgrid(r,theta)

#Cloud Ellipse parameters
b = rc_par # semi-major axis
bc = 0.6 #approximate braginskii coefficient 
a = bc*b #semi-minor axis

#Pellet Ellips parameters

bp = rp_par
bcp = 0.6
ap = bp/bcp

all_rc = np.zeros(len(theta))
all_rp = np.zeros(len(theta))
rc_ind = np.zeros(len(theta))

all_rp[:] = rp_par # if you want a spherical pellet and elliptical cloud just use this

D = np.zeros(np.shape(R))
for i in range(0, len(theta)):
    all_rc[i] = a*b/(np.sqrt((a*np.cos(theta[i]))**2 + (b*np.sin(theta[i]))**2))
    #all_rp[:] = ap*bp/(np.sqrt(ap*np.cos(theta[i])**2) + (bp*np.sin(theta[i])**2))
    rc_ind[i] = (np.abs(r[:] - all_rc[i])).argmin()
    D[i,ind_ppar:int(rc_ind[i])] = 0.01*((1.0 + all_rp[i]**2)/(1.0 + R[i,ind_ppar:int(rc_ind[i])])**2)

cmap = plt.get_cmap('viridis')
ax.set(xlim = (0.0,2*np.pi), ylim = (0.0,35.0))
ax.set_yticks([5,15,25,35])
ax.set_xticks([])
#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
im = ax.pcolormesh(T, R, D, cmap=cmap, norm = colors.LogNorm(vmax = 0.01, vmin = D[0,int(rc_ind[0])-1]))
plt.grid(axis = 'y', color = 'black')
fig.colorbar(im, ax=ax)
plt.text(45,40, r'$\rightarrow \mathrm{\mathbf{B}}$', fontsize = 14)
plt.savefig('elliptical_cloud.png', format = 'png', dpi = 1200)
plt.show()