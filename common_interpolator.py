import numpy as np 
import scipy.interpolate as spint 

"""This function adds A2 defined on the points r2 to A1 defined on the points r1
Must make sure they are defined on all common points on the common grid 
to be added - this will be relevant for the electron densities that have stopped
and the electrons that have been transported also.
"""
def common_interp(r, A1, r1): #make sure r1 and r2 are in the common grid with pellet centre at zero of the coord system 

    ind_interp_low = next(p[0] for p in enumerate(r) if p[1] >= r1[-1]) # lowest common point
    ind_interp_up = next(p[0] for p in enumerate(r) if p[1] > r1[0]) # greatest common point

    g = spint.interp1d(r1, A1, kind = 'cubic') # defining interpolation function

    A1_new = g(r[ind_interp_low:ind_interp_up]) # interpolation across the defined range
    return A1_new, ind_interp_low, ind_interp_up

    