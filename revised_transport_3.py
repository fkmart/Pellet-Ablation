import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rc, rp, eps, dt, t 
from stop_calc_rp_rc import rdot, rcdot 

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
for i in range(0, len(trans_r)):
    if trans_r[i] > rc2:
        arr = np.append(arr,i)

trans_r = np.delete(trans_r, arr)
trans_dens = np.delete(dens1, arr)

"Now define dens2 on trans_r points"

dens1_check = eps*(1.0 + rp1**2)/(1.0 + trans_r**2)
dens2_check = eps*(1.0 + rp2**2)/(1.0 + trans_r**2)

fracs = dens2_check/trans_dens

fig, ax = plt.subplots()
ax.plot(r1, dens1)
ax.plot(r2, dens2)
ax.plot(trans_r,(1.0 - fracs)*dens1_check)
ax.plot(trans_r, fracs*trans_dens)

"Now add the two parts of the outer profile together"
#ax.plot(trans_r, (1.0 -fracs)*dens1_check + fracs*trans_dens)
plt.show()