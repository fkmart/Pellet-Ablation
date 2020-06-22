import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp, rc, rpd, t , eps
from u_cloud_sol import u_sol 
import romberg as romb
import scipy.interpolate as spint
"""
#SETUP VARIABLES
t1 = 40
t2 = 50 

shift = t2 -t1 

rp1, rp2 = rp[t1], rp[t2]
rc1, rc2 = rc[t1], rc[t2] 

l = 1025

r = np.linspace(rp1, rc2, num = l, endpoint = 'true')

ind_rc1 = next(p[0] for p in enumerate(r) if p[1] > rc1)

rho1 = eps*(1.0 + rp1**2)/(1.0 + r[:ind_rc1]**2)

#MASS TEST VARIABLES 
r1 = np.linspace(rp1, rc1, num = l, endpoint = 'true')
dens1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)

u = u_sol(r[:ind_rc1], t1)

f_rem = (1.0 + rp2**2)/(1.0 + rp1**2)
f_tran = 1.0 - f_rem 

r_trans = r[:ind_rc1] + u*shift*(t[1] - t[0])

ind_low_trans = next(p[0] for p in enumerate(r) if p[1] > r_trans[0])

#INTERPOLATION
f = spint.interp1d(r_trans, f_tran*rho1, fill_value = 'extrapolate')
rho_trans_interp = f(r[ind_low_trans:])

dens_final = np.zeros(l)
dens_final[:ind_rc1] += f_rem*rho1
dens_final[ind_low_trans:] +=  rho_trans_interp 

#MASS CONSERVATION CHECK 
print('before transport : ' + str(romb.romberg_samp(dens1,r1)))
print('after transport : ' + str(romb.romberg_samp(dens_final, r)))"""

def transport(t1,shift,rgrid, density):
    from gen_var import rp,rc 
    t2 = t1 + shift
    rp1,rc1 = rp[t1], rc[t1]
    rp2,rc2 = rp[t2], rc[t2]
    
    u = u_sol(rgrid, t1)
    f_rem = (1.0 + rp2**2)/(1.0 + rp1**2)
    f_tran = 1.0 - f_rem 

    r_trans = r + u*shift*(t[1] - t[0])

    dens_rem = np.asarray([density*f_rem,r])
    dens_trans = np.asarray([density*f_trans, r_trans])
    return dens_rem, dens_trans
