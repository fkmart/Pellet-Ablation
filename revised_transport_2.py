import numpy as np 
import matplotlib.pyplot as plt 
import math as mt 
from gen_var import rc,rp, eps, t , dt
from stop_calc_rp_rc import rdot, rcdot 
import grid_pusher_3 as gp 
import array_delete as ad


t1 = 40 
t2 = 60 
shift = t2-t1 

rp1 = rp[t1] 
rp2 = rp[t2]
rc1 = rc[t1] 
rc2 = rc[t2] 

l = 2000

r1 = np.linspace(rp1,rc1, num = l, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = l, endpoint = 'true')

grid = np.linspace(rp2,rc2,num = 10000, endpoint = 'true')

dens1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)
dens2 = eps*(1.0 + rp2**2)/(1.0 + r2**2)

frac = dens2[-1]/dens1[-1]
print(frac)

"Now calculate the fluid velocities"

rdots = rdot(t,rc,rp)
rcd = rdots[t1]

const = dens1[-1]*rcd

u1 = const/dens1[:]

trans_r = r1 + dt*shift*u1
trans_dens1 = frac*dens1

dens_grid = gp.pusher(trans_r, trans_dens1, r2)
pushed_fdens1 = gp.pusher(r1, (1-frac)*dens1, r2)

dens_grid += pushed_fdens1
#dens_grid += pushed_fdens1

arr = ad.array_delete_less(dens_grid, 1e-14)

dens_grid = np.delete(dens_grid, arr)
grid = np.delete(r2, arr)


#arr = ad.array_delete_more(dens_grid)

"Plotting"
fig, ax  = plt.subplots() 
ax.plot(r1,dens1, label = r'$\tilde{\rho}_1$')
ax.plot(r2, dens2, label = r'$\tilde{\rho}_2$')
ax.plot(trans_r, trans_dens1, label = r'$\tilde{\rho}_{1 \rightarrow 2}$')
ax.plot(r1, (1-frac)*dens1, label = r'$(1-f) \tilde{\rho}_1$')
#ax.plot(r2, dens_grid, label = 'new dens')
#ax.plot(grid_final, dens_final, label = r'$\tilde{\rho}_{\mathrm{new}}$')
plt.legend()
plt.show()