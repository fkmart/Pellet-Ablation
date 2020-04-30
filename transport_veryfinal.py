import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rc, rp, t,dt  , eps 
from stop_calc_rp_rc import rdot, rcdot 
from transport_functions import transport
import scipy.interpolate as spint

t1 = 40 
t2 = 55 
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

"Now calculate the fluid velocities"

rdots = rdot(t,rc,rp)
rcd = rdots[t1]

const = dens1[-1]*rcd

u1 = const/dens1[:]
trans_r = r1 + dt*shift*u1

"Now iterate through the function to get all fractional mass motions"

k = 10.0

fig,ax = plt.subplots() 
dexp = eps*np.exp(-k*(r2 -rp2)**2)
ax.plot(r2,dexp)
plt.yscale('log')
plt.show()
ind_low = next(p[0] for p in enumerate(r2) if p[1] > rc1) + 1
r_test = r2[ind_low]
ind_up = l
all_frac_r = [] 
all_frac_t = [] 
frac_rem = np.zeros(np.shape(r2[ind_low:ind_up]))
r_plot = r2[ind_low:ind_up]
dens_new = eps*(1.0 + rp2**2)/(1.0 + r_plot**2)
dens_exp = eps*np.exp(-k*(r_plot - rp2)**1)
frac_trans, pre_r = transport(r1,r2,trans_r, k,frac_rem, ind_low, ind_up)
dens_old = eps*(1.0 + rp1**2)/(1.0 + pre_r**2)
frac_rem = 1.0 - frac_trans

#Plotting block was tested - returns correct result for new cloud area
"""fig, ax = plt.subplots() 
ax.plot(r_plot, dens_new)
ax.plot(r_plot, dens_exp)
ax.plot(r_plot, frac_rem*dens_old)
ax.plot(r_plot, dens_exp + frac_trans*dens_old, linestyle = '--')
plt.show()"""
counter = 0
while pre_r[0] > trans_r[0]:
    if counter ==14:
        break
    else:
        pass
    counter +=1
    all_frac_r.append(np.flip(frac_rem, axis = 0))
    all_frac_t.append(np.flip(frac_trans, axis = 0))
    ind_up = np.copy(ind_low)
    ind_low = next(p[0] for p in enumerate(r2) if p[1] >pre_r[0])
    h = spint.interp1d(pre_r, frac_rem, kind = 'cubic', fill_value = 'extrapolate') #interpolate fractions to proper points
    frac_rem_int = h(r2[ind_low:ind_up])
    frac_trans, pre_r = transport(r1,r2,trans_r,k,frac_rem_int, ind_low,ind_up)
    frac_rem = 1.0 - frac_trans
    #define densities to check
    r_plot = r2[ind_low:ind_up]
    dens_new = eps*(1.0 + rp2**2)/(1.0 + r_plot**2)
    dens_exp = eps*np.exp(-k*(r_plot - rp2)**2)
    dens_trans = eps*(1.0 + rp1**2)/(1.0 + pre_r**2)
    dens_old = eps*(1.0 + rp1**2)/(1.0 + r_plot**2)
    fig, ax  = plt.subplots()
    ax.plot(r_plot, dens_new, label = 'new dens')
    ax.plot(r_plot, dens_old, label = 'old dens')
    ax.plot(r_plot, dens_exp, label = 'exponential')
    ax.plot(r_plot, frac_trans*dens_trans, label = 'transported')
    ax.plot(r_plot, frac_rem_int*dens_old, label = 'remainder')
    ax.plot(r_plot, dens_exp + frac_rem_int*dens_old + frac_trans*dens_trans, linestyle = '--', label = 'total')
    plt.legend()
    plt.show()
