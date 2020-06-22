import numpy as np 
import scipy.interpolate as spint 

def u_source(t,rp,rc,r, alpha):
    from gen_var import rp,rc, eps 
    from stop_calc_rp_rc import rdot
    import math as mt 
    rp1 = rp[t]
    rc1 = rc[t]
    dens = eps*(1.0 + rp1**2)/(1.0 + r**2)

    rpd = -0.5*((np.log(np.cosh(1.0)))**(-0.5))*(((np.log(np.cosh(1.0 - t))))**(-0.5))*(np.tanh(1-t))
    rd1 = rdot(t,rc1,rp1)
    rcd = rd1[t,-1]
    #Now to input the analytical solution of the source/sink continuity equation
    
    u = (alpha + dens[-1]*rcd - 2.0*eps*rpd*rp*(mt.arctan(r) - mt.arctan(rc1)))/(dens)
    return u 