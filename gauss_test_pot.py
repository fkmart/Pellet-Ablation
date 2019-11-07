import numpy as np 
import scipy.interpolate as spint

def gauss_func(A,k,x0,x):
    pot = A*np.exp(k*(x - x0)**2)
    return pot 

def faux_pot(r,dens,A, r_full):
    f = spint.interp1d(r,dens, kind = 'cubic')
    ind_low = next(p[0] for p in enumerate(r_full) if p[1] > r[0])
    ind_up = next(p[0] for p in enumerate(r_full) if p[1] > r[-1])
    faux_pot = f(r_full[ind_low:ind_up])
    faux_pot *= A # scaling factor to get a reasonable amount of slowing 
    return faux_pot 

