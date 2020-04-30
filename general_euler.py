import numpy as np 
import matplotlib.pyplot as plt 
import rkf45_euler 
import scipy.interpolate as spinter
from gen_var import rp,rc, eps, dt, t
import math as mt 
from stop_calc_rp_rc import rcdot, rdot 

t1 = 40
t2 = 50 

rp1 = rp[t1]
rp2 = rp[t2]

rc1 = rc[t1] 
rc2 = rc[t2]

r1 = np.linspace(rp1,rc1, num = 400, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = 400, endpoint = 'true')

dens1 = eps*((1.0 + rp1**2)/(1.0 + r1[:]**2))
dens2 = eps*((1.0 + rp2**2)/(1.0 + r2[:]**2)) 

rpd1 = 0.5*mt.tanh(1.0 - t1*dt)/(np.sqrt(np.log(mt.cosh(1.0)))*np.sqrt(np.log(mt.cosh(1 - t1*dt))))
def euler(r,rd):
    rdr = 2.0*r*rd/(1 + r**2) - 2.0*rpd1*rp1/(1 + rp1**2)
    return rdr

#butcher-tableau 
BT = np.zeros((6,6))
BT[1,:] = [0.25,0.25,0.0, 0.0, 0.0, 0.0]
BT[2,:] = [0.125, 3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0]
BT[3,:] = [12.0/13.0, 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0]
BT[4,:] = [1.0, 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0]
BT[5,:] = [0.5, -8.0/27.0, 2.0, -3544.0/2565.0, 1859./4140.0 , -11.0/40.0]
    
#butcher-tableau for solutions

BTS = np.zeros((2,6))
BTS[0,:] = [25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, - 0.2, 0.0]
BTS[1,:] = [16.0/135.0, 0.0, 6656.0/12825.0,28561.0/56430.0, -9.0/50.0, 2.0/55.0 ]

rcdots= rdot(t,rc,rp) 
rcd1 = rcdots[t1]
h = 1e-2

rdots, not_needed = rcdot(t,rc,rp)

rd_plot = rdots[t1,:]
r_plot = np.linspace(rp1,rc1, len(rd_plot), endpoint = 'true')
r = np.asarray([rc1])
rd = np.asarray([rcd1])
h_arr = np.asarray([h])
err_arr = []

while r[-1] > rp1:
    r_new, rd_new, h,err,ks  = rkf45_euler.rkf(r[-1],rd[-1],h,euler,BT,BTS, rc1,rp1)
    r = np.append(r,r_new)
    rd = np.append(rd,rd_new)
    h_arr = np.append(h_arr,h)
    err_arr = np.append(err_arr, err)
print('done!')
fig, ax = plt.subplots() 
ax.plot(r,rd)
ax.plot(r_plot, rd_plot)
plt.show()