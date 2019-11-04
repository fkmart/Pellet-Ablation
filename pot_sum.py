import numpy as np

def pot_sum(elec_pot, pot,elec_r, pot_r):
    low_ind = np.argmin(np.abs(elec_pot[:] - pot_r[0])) # lowest matching index
    up_ind = np.argmin(np.abs(elec_pot[:] - elec_r[-1])) #highest matching index
    tot_pot = elec_pot + pot[low_ind:up_ind] # add potentials at corresponding r-values
    return tot_pot