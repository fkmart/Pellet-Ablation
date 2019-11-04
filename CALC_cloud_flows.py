import numpy as np
import matplotlib.pyplot as plt
import stop_calc_rp_rc
from gen_var import dr, eps, t, dt
import math as mt 
import scipy.interpolate as spint

rp = stop_calc_rp_rc.calc_rp(t)
rc = stop_calc_rp_rc.calc_rc(t)
rcd = stop_calc_rp_rc.rdot(t,rc,rp)
 
t1 = 450
t2 = 451
delta_t = t2 - t1

dr = 0.5
#IMPORTANT COMMENT 
#Can not simply start FROM RP and step in dr to RC to get an array for r
#Should instead START AT RC and move BACKWARDS by dr to RP
"for t1"
rp1 = rp[t1]
rc1 = rc[t1]
rcd1 = rcd[t1]
r1 = np.arange(rc1,rp1, -dr)
dens1 = eps*((1.0 + rp1**2)/(1.0 + r1**2))

"for t2"
rp2 = rp[t2] 
rc2 = rc[t2]
rcd2 = rcd[t2] 
r2 = np.arange(rc2, rp2, -dr)
dens2 = eps*((1.0 + rp2**2)/(1.0 + r2**2))

new_rc = rcd1*(delta_t*dt) + rc1
difference = rc2 - new_rc
print('DIfference in r_c co-ordinates is ' + str(difference))

#Need to define a new grid that exists on itneger multiples of dr so that
#one grid can easily be interpolated onto the other to have some common frame for the 
#whole system. For now, let's just evaluate the internal points and verify that 
#the earlier time comfortably advances onto the next cloud

#THIS IS KIND OF RIGHT BUT LET'S TRY SOMETHING DIFFERENT INSTEAD
"New grid for t1"
"""
rpn1 = rp1 + dr   
rcn1 = rc1 - dr 
indp = mt.floor(rpn1/dr)
indc = mt.floor(rcn1/dr)

r1_newgrid = np.arange(indp*dr, indc*dr + dr, dr)
f1 = spint.interp1d(r1,dens1, kind = 'cubic')
dens1_new = f1(r1_newgrid)

"NEw grid for t2"

rpn2 = rp2 + dr
rcn2 = rc2 - dr    
indp2 = mt.floor(rpn2/dr)
indc2 = mt.floor(rcn2/dr) 

r2_newgrid = np.arange(indp2*dr, indc2*dr + dr, dr) 
f2 = spint.interp1d(r2, dens2, kind = 'cubic')
dens2_new = f2(r2_newgrid)
"""
###########################################################
###########################################################

"""Instead do the following:
1) Calculate internal speeds on regular grid
2) Calculate factor difference in rc1 to rc2 in density
3)Apply everywhere else and see what happens"""

r1_int = dens1[0]*rcd1/dens1[:]
r1_2 = r1[:] + (dt*delta_t*r1_int[:])
dens1_2 = dens1[:]


fig, ax = plt.subplots()

ax.plot(r1, dens1, marker = 'x' ,label = '1',linestyle = '')
ax.plot(r2, dens2, label = '2', marker = 'x', linestyle = '')
ax.plot(r1_2, dens1_2, label = '1 - 2 old',marker =  'x', linestyle = '')
#ax.plot(r2,dens2)
#ax.plot(r1, dens1)
#ax.plot(r2_newgrid, dens2_new)
#ax.plot(r1_newgrid, dens1_new)
plt.legend()
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel(r'$\tilde{\rho}}$', rotation = 0)
ax.set_xlabel(r'$\tilde{r}$')
plt.show()
print('something celebratory')