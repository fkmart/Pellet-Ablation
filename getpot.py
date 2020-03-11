import numpy as np 
import scipy.interpolate as spit 

def getpot(pot,h,x_start):
    from gen_var import r
    x_new = x_start - h 
    ind_start = (r[:] - x_start).argmin()
    ind_end = (r[:] - x_new).argmin()
    pot_diff = pot[ind_start] - pot[ind_end]
    return pot_diff