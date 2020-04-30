import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rc, rp, eps, dt, t 
from stop_calc_rp_rc import rdot, rcdot 
import scipy.interpolate as spint 

t1 = 40 
t2 = 45 
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

arr = []
trans_r = r1 + dt*shift*u1
"""for i in range(0, len(trans_r)):
    if trans_r[i] > rc2:
        arr = np.append(arr,i)

trans_r = np.delete(trans_r, arr)
trans_dens = np.delete(dens1, arr)"""

"Now define dens2 on trans_r points"

dens1_check = eps*(1.0 + rp1**2)/(1.0 + trans_r**2)
dens2_check = eps*(1.0 + rp2**2)/(1.0 + trans_r**2)

"Exponential Input"

k = 1.0 # subject to change 
dens_exp = 0.01*np.exp(-k*(r2 - rp2))

"Now need to interpolate trans_r points to r2 points"

r_range = np.linspace(trans_r[0], rc2, num = len(trans_r), endpoint = 'true')
f = spint.interp1d(trans_r, r1, kind = 'cubic')
pre_transit_r = f(r_range)

pre_trans_dens = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)
d = eps*(1.0 + rp2**2)/(1.0 + r_range**2)
d_exp = 0.01*np.exp(-k*(r_range - rp2))
f_trans = (d - d_exp)/pre_trans_dens

fig,ax = plt.subplots() 
ax.plot(r_range, d)
ax.plot(r_range, d_exp)
ax.plot(r_range, f_trans*pre_trans_dens)
ax.plot(r_range, f_trans*pre_trans_dens + d_exp)
plt.show()
print(f_trans)
f_rem = 1.0 - f_trans 

#Things seem fine on an initial look
#Now need to interpolate to the same points

ind1 = next(p[0] for p in enumerate(r2) if p[1] > pre_transit_r[0])
ind2 = next(p[0] for p in enumerate(r2) if p[1] > r_range[0])

r_new1 = r2[ind1:]
r_new2 = r2[ind2:]

g = spint.interp1d(pre_transit_r, f_rem*dens1, kind = 'cubic', fill_value = 'extrapolate')
dens_new1 = g(r_new1)

h = spint.interp1d(r_range, f_trans*pre_trans_dens, kind = 'cubic', fill_value = 'extrapolate')
dens_new2 = h(r_new2)

dens_final = np.copy(dens_exp)
dens_final[ind1:] += dens_new1
dens_final[ind2:] += dens_new2
fig,ax = plt.subplots()
ax.plot(r2, dens_exp, label = 'exponential input')
ax.plot(pre_transit_r, f_rem*dens1, label = 'remaining dens')
ax.plot(r_range,f_trans*pre_trans_dens,label= 'transported dens')
ax.plot(r2,dens2, label = 'new dens', linestyle = '--')

ax.plot(r2, dens_final, label = 'final density')
plt.legend()
plt.show()