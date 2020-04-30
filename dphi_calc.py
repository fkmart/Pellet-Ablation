import numpy as np
import scipy.interpolate as spint 

def dphi_calc(x1,x2,r, pot):
    int_f = spint.interp1d(r, pot, kind = 'cubic', fill_value = 'extrapolate')
    last_pot = pot[-1]
    rc_thing = r[-1]
    pot1 = int_f(np.asarray(x1))
    pot2 = int_f(np.asarray(x2))
    dphi = pot2 - pot1 
    return dphi