import numpy as np 
import scipy.integrate as spint 
from gen_var import dr,pel_pot, rp, rc, t_start, r
import scipy.interpolate as spit
import os 
import romberg as ro 

def renorm(r_raw, bins_raw):
    frac_left = np.sum(bins_raw)
    #Interpolation Information
    rpi = rp[t_start]
    rci = rc[t_start]
    f = spit.interp1d(np.flip(r_raw,axis = 0), np.flip(bins_raw, axis = 0), 
        fill_value = 'extrapolate', kind = 'quadratic')
    g = spit.PchipInterpolator(np.flip(r_raw, axis = 0), np.flip(bins_raw, axis = 0), extrapolate = 'true')
    low_int = next(p[0] for p in enumerate(r) if p[1] > r_raw[-1]) - 1
    up_int = next(p[0] for p in enumerate(r) if p[1] > r_raw[0])

    low = next(p[0] for p in enumerate(r) if p[1] > rpi)
    up = next(p[0] for p in enumerate(r) if p[1] > rci)

    r1 = r[low_int]
    r2 = r[up_int]
    r_int = np.linspace(rpi,rci, num = 1 + 2**10, endpoint = 'True')
    n_romb = len(r)
    bins_full = f(r_int)
    integrated = spint.romb(bins_full, dx = -r_int[0] + r_int[1])
    dens_full = bins_full[:]/integrated
    return bins_full, dens_full, r_int