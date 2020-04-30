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
ind_cloud = next(p[0] for p in enumerate(r2) if p[1] > rc1)

r_range = np.linspace(r2[ind_cloud], rc2, num = len(trans_r), endpoint = 'true')
f = spint.interp1d(trans_r, r1, kind = 'cubic')
pre_transit_r = f(r_range)

pre_trans_dens = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)
d = eps*(1.0 + rp2**2)/(1.0 + r_range**2)
d_exp = 0.01*np.exp(-k*(r_range - rp2))
f_trans = (d - d_exp)/pre_trans_dens
f_rem = 1.0 - f_trans 

all_f_trans = []
all_f_rem = [] 
all_f_trans.append(f_trans)
all_f_rem.append(f_rem)
"""
fig,ax = plt.subplots() 
ax.plot(r_range, d)
ax.plot(r_range, d_exp)
ax.plot(r_range, f_trans*pre_trans_dens)
ax.plot(r_range, f_trans*pre_trans_dens + d_exp, linestyle = '--')
plt.show()
""" 

"""Now need to shift this function back into common points and
also consider the stuff that moves INTO those points"""
dens2 = eps*(1.0 + rp2**2)/(1.0 + pre_transit_r**2)
dens1 = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)

ind_int = next(p[0] for p in enumerate(r2) if p[1] >pre_transit_r[0])
g = spint.interp1d(r2,pre_transit_r, kind = 'cubic')
pre_transit_r2 = g(pre_transit_r)
pre_trans_dens2 = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r2**2)
d2 = eps*(1.0 + rp2**2)/(1.0 + pre_transit_r**2)
d1 = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)
d_exp = eps*np.exp(-k*(pre_transit_r - rp2))
f_trans_new = (d2 - d_exp - f_rem*d1)/pre_trans_dens2 

fig,ax = plt.subplots()
ax.plot(pre_transit_r, d2)
ax.plot(pre_transit_r, d_exp)
ax.plot(pre_transit_r, f_trans_new*pre_trans_dens2)
ax.plot(pre_transit_r, f_rem*d1)
ax.plot(pre_transit_r, d_exp + f_rem*d1 + f_trans_new*pre_trans_dens2, linestyle = '--')
plt.show()

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