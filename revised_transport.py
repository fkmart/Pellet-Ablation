import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as spint 
from gen_var import rc, rp , eps, dt
from stop_calc_rp_rc import rdot 
from gen_var import t as time
import grid_pusher_3 as gp 
import scipy.interpolate as spint


t1 = 40
t2 = 50 

shift = t2-t1

rp1 = rp[t1]
rp2 = rp[t2] 

rc1 = rc[t1] 
rc2 = rc[t2] 

r_input = np.linspace(rp2,rc2,num = 500, endpoint = 'true')

#equilibrium profiles
r_eq_1 = np.linspace(rp1,rc1, num = 500, endpoint = 'true')
r_eq_2 = np.linspace(rp2,rc2, num = 500, endpoint = 'true')

dens_eq_1 = eps*((1.0 + rp1**2)/(1.0 + r_eq_1**2))
dens_eq_2 = eps*((1.0 + rp2**2)/(1.0 + r_eq_2**2))

#input dens from ablation 
k = -np.log((rp1**3 - rp2**3)/(3.0*eps*(rp2*(1.0 - rp2) + rc2*(rc2 - 1.0))))
dens_in = eps*np.exp(-k*(r_input - rp2))
#Now look at equilibrium value at rc2, subtract the input exponential

rc2_rem = dens_eq_2[-1] - dens_in[-1]

frac = rc2_rem/dens_eq_1[-1]

print(frac)

#Now transport this fraction forward

rdots = rdot(time,rc,rp)
rcd = rdots[t1]

const = dens_eq_1[-1]*rcd

rp1_ind = next(p[0] for p in enumerate(r_eq_1) if p[1] > rp1)

transport_mass = np.zeros(len(r_eq_1))
transport_mass[:] = dens_eq_1[:]*frac
rd = const/dens_eq_1[:]
trans_r = r_eq_1[:] + dt*rd[:]*shift

#Now need to define a romberg grid to push densities to 

grid = np.linspace(rp2,rc2, endpoint = 'true', num = 2049)
all_dens = np.append(dens_in,(1-frac)*dens_eq_1)
all_dens = np.append(all_dens, frac*dens_eq_1)

all_r = np.append(r_input, r_eq_1)
all_r = np.append(all_r, trans_r)

push_dens = gp.pusher(all_r,all_dens, grid)

#deletes zeroes in the pushed array

arr = []
for i in range(0, len(push_dens)):
    if push_dens[i] < 1e-14:
        arr = np.append(arr,i)
    else:
        pass
old_dens = np.delete(push_dens, arr)
old_grid = np.delete(grid, arr)

g = spint.interp1d(old_grid, old_dens, kind = 'cubic')
new_dens = spint.pchip_interpolate(old_grid, old_dens, grid)
new_dens = g(grid)
fig,ax = plt.subplots()
#ax.plot(r_input,dens_in)
ax.plot(trans_r,transport_mass)
#ax.plot(r_eq_1,(1-frac)*dens_eq_1)
ax.plot(r_eq_2,dens_eq_2)
#ax.plot(grid, push_dens)
plt.show()