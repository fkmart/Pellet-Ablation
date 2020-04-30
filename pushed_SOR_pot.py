import numpy as np
import matplotlib.pyplot as plt 
import grid_pusher as gp 
import os  
from gen_var import rp,rc 
import iterative_sol as SOR
import discret_mat as dm
import scipy.interpolate as spint
import romberg as ro

"""This code loads in a CSDA and a CSDA-sheath run to compare distribution of charge 
and calculate the potential. 

Procedure is the same for both:
1) Densities are interpolated and a romberg integrator determines enclosed area

2) Zeroes on the grid are removed and remaining values are pushed to a grid.

3)Results are interpolated again onto the grid post-push and integrated via romberg again

4)Densities are renormalised such that enclosed area yields same value as pre-push

5)SOR Solver determines the potential on this grid."""

direc = os.getcwd()
t = 40 
load_dir = os.path.join(direc, 'one_iteration', 'analysed_outputs', 't_' + str(t)) + os.sep 

file = np.loadtxt(load_dir + 'density_neutral.txt')
dens = file[:,0]
pos = file[:,1]

file2 = np.loadtxt(load_dir + 'density_sheath.txt')
dens2 = file2[:,0]
pos2 = file2[:,1]

pos = np.flip(pos, axis = 0)
dens = np.flip(dens, axis = 0)

pos2 = np.flip(pos2, axis = 0)
dens2 = np.flip(dens2, axis = 0)
res = 257
r_grid = np.linspace(rp[t],rc[t], res, endpoint = 'true')
print('Difference between points is ' + str(r_grid[1] - r_grid[0]))

def func(dens,pos, r_grid,fi):
    dens_interpolated = spint.pchip_interpolate(pos,dens,r_grid)

    integrated_dens = ro.romberg_samp(dens_interpolated, r_grid)

    push_dens = np.zeros(res)
    push_dens = gp.pusher(fi, r_grid)

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
    brand_new_dens = spint.pchip_interpolate(old_grid,new_dens,r_grid)
    integrated_bnd = ro.romberg_samp(brand_new_dens,r_grid)
    norm = integrated_dens/integrated_bnd
    brand_new_dens *= norm
    print('Renormalisation constant is '+ str(norm))
    print(integrated_dens)
    print(ro.romberg_samp(brand_new_dens,r_grid))
    return brand_new_dens, new_dens, old_grid, push_dens

bnd_neut, pushed_neut, og_neut, pdn = func(dens,pos,r_grid,file)
bnd_char, pushed_char, og_char, pdc = func(dens2,pos2,r_grid,file2)
#brand_new_dens = f(r_grid)
#s = np.sum(brand_new_dens)
#brand_new_dens /=s

A = dm.discret(r_grid)

x0 = 1e-3
e0 = 1.6*10**(-19)
n0 = 1e19
eps0 = 8.854*10**(-12)
V0 = 1e3
const = x0*x0*e0*n0/(eps0*V0)

pot = np.zeros(res)
pot[0] = 0.0
#pot = SOR.SOR(A,pot,const*push_dens,r_grid)
pot_neut = SOR.SOR(A,pot,bnd_neut*const,r_grid)

pot_ch = np.linspace(-2.8,0.0, num = len(r_grid), endpoint = 'true')
pot_ch = SOR.SOR(A,pot_ch,bnd_neut*const,r_grid)

fig,ax = plt.subplots(figsize = (8,6))
plt.text(1.0,-2.0 ,r'$\leftarrow$'  + 'To pellet')
l1 = ax.plot(r_grid, pot, color = 'navy', label = r'$\tilde{\phi}$')
ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12, rotation = 0)
ax.set_xlabel(r'$\tilde{r}$')
ax2 = ax.twinx() 
l2 = ax2.plot(r_grid, bnd_neut, color = 'navy',linestyle = '--', label = r'$\tilde{n}$')
ax2.set_ylabel(r'$\tilde{n}$', fontsize = 12, rotation = 0)

l3 = ax.plot(r_grid, pot_ch, color = 'orange', label = r'$\tilde{\phi}_c$')
l4 = ax2.plot(r_grid, bnd_char, color = 'orange', linestyle = '--', label = r'$\tilde{n}_c$')

lines = l1 + l2 +l3 + l4
labs = [l.get_label() for l in lines]
plt.legend(lines, labs, loc = 'lower right')
plt.title('1D Poisson solution on charged and neutral pellets')

plt.tight_layout()
plt.savefig('SOR_sheathed_pellet.png', format='png', dpi = 1400)
plt.show()
"""
fig,ax = plt.subplots()
ax.plot(pos,dens, color = 'blue', label = 'Raw + neutral')
ax.plot(pos2,dens2, color = 'orange', label = 'Raw + charged')
ax.plot(r_grid, bnd_neut, color = 'green', label = 'Pushed + neutral')
ax.plot(r_grid, bnd_char, color = 'purple', label = 'Pushed + charged')
plt.legend()
plt.show()"""