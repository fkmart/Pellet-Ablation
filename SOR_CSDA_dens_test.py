import numpy as np 
import os 
import matplotlib.pyplot as plt 
import scipy.interpolate as spint 
import scipy.integrate as spit
from gen_var import rp,rc,r 
import iterative_sol as SOR
import discret_mat

cd = os.getcwd() 
path = os.path.join(cd, 'one_iteration','analysed_outputs')
t = 50 
load_path = os.path.join(path, 't_' + str(t))
file = np.loadtxt(load_path + os.sep + 'density.txt')
dens = file[:,1]
pos = file[:,0] 
low = next(p[0] for p in enumerate(r) if p[1] > rp[t]) 
up = next(p[0] for p in enumerate(r) if p[1] > rc[t]) 
r_cloud = r[low:up]

g = spint.interp1d(pos,dens, kind = 'quadratic', fill_value = 'extrapolate')
dens_full = g(r_cloud)
sum_dens = np.sum(dens) 
renorm_dens = sum_dens*dens_full/(np.sum(dens_full))

x0 = 1e-3
q0 = 1.6*10**(-19)
n0 = 10**19
esp0 = 8.854*10**(-12)
V0 = 1e3
K = x0*x0*q0*n0/(esp0*V0)
print(K)

A = discret_mat.discret(r_cloud)
phi = np.zeros(len(r_cloud))
phi[0] = 0.0
phi[-1] = 0.0
phi = SOR.SOR(A,phi, K*dens_full,r_cloud)

fig, ax  = plt.subplots() 
ax2 = ax.twinx() 
l1 = ax.plot(r_cloud,phi*V0, color = 'blue', label = r'$\tilde{\phi}$')
l2 = ax2.plot(r_cloud, dens_full, color = 'red', label = r'$\tilde{n}$')
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0)
ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12, rotation = 0)
ax2.set_ylabel(r'$\tilde{n}$', fontsize = 12, rotation = 0)
lines = l1 + l2 
labs = [l.get_label() for l in lines]
ax.legend(lines,labs)
#plt.savefig('SOR_CSDA_test.png', format = 'png', dpi = 1200)
plt.show()