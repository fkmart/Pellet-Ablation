import numpy as np
from stop_calc_rp_rc import *
import matplotlib.pyplot as plt

rp = calc_rp(t)
print('rp done')
rc = calc_rc(t)
print('rc done')
rcd, r_i = rdot(t,rc,rp)
print('rcd done')

t1 = 100
t2 = 120
dt = t2 - t1 
lr = 100

r1 = np.linspace(rp[t1], rc[t1], lr)
r2 = np.linspace(rp[t2], rc[t2], lr)

rcd1 = rcd[t1]
rcd2 = rcd[t2]

def density(r):
    dens = eps*((1.0 + r[0]**2)/(1.0 + r[:]**2))
    return dens 

d1 = density(r1)
d2 = density(r2) 

cons = 0.5*(d1[-1]*rcd1[0] + d2[-1]*rcd2[0])

rd = np.zeros(lr)

for i in range(0, lr):
    rd[i] = cons/d1[i]