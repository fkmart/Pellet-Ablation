import numpy as np
import scipy.interpolate as spint 
import makeshift_interpolator as MI
from numba import njit

def dphi_calc(x1,x2,r, pot):
    int_f = spint.interp1d(r, pot, kind = 'cubic', fill_value = 'extrapolate')
    pot1 = int_f(np.asarray(x1))
    pot2 = int_f(np.asarray(x2))
    dphi = pot2 - pot1 
    return dphi

@njit
def dphi_calc_jit(x1,x2,r, pot):
    pot1 = MI.makeshift(x1,r,pot)
    pot2 = MI.makeshift(x2,r,pot)
    dphi = pot2 - pot1 
    return dphi