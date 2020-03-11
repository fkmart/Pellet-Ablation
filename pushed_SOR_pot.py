import numpy as np
import matplotlib.pyplot as plt 
import grid_pusher as gp 
import os  
from gen_var import rp,rc 
import iterative_sol as SOR
import discret_mat as dm
import scipy.interpolate as spint

direc = os.getcwd()
t = 40 
load_dir = os.path.join(direc, 'one_iteration', 'analysed_outputs', 't_' + str(t)) + os.sep 

file = np.loadtxt(load_dir + 'density.txt')
dens = file[:,0]
pos = file[:,1]

res = 256
r_grid = np.linspace(rp[t],rc[t], res, endpoint = 'true')
push_dens = np.zeros(res)
push_dens = gp.pusher(file, r_grid)
"""
#Following block of code removes zero elements from pushed dens and interpolates
arr = []
for i in range(0, res):
    if push_dens[i] < 1e-14:
        arr = np.append(arr,i)
    else:
        pass


new_dens = np.delete(push_dens, arr)
old_grid = np.delete(r_grid, arr)

f = spint.interp1d(old_grid,new_dens,kind = 'cubic', fill_value = 'extrapolate')
brand_new_dens = f(r_grid)
s = np.sum(brand_new_dens)
brand_new_dens /=s

"""
A = dm.discret(r_grid)

x0 = 1e-3
e0 = 1.6*10**(-19)
n0 = 1e19
eps0 = 8.854*10**(-12)
V0 = 1e3
const = x0*x0*e0*n0/(eps0*V0)

pot = np.zeros(res)
pot[0] = -2.8
pot = SOR.SOR(A,pot,const*push_dens,r_grid)
#pot = SOR.SOR(A,pot,2.0*brand_new_dens*const,r_grid)

fig,ax = plt.subplots() 
ax.plot(r_grid, pot, color = 'navy')
ax2 = ax.twinx() 
ax2.plot(r_grid, push_dens, color = 'orange')
plt.legend()
plt.show()