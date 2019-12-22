import numpy as np
import os
from gen_var import *
from electron import dist_calc, e_bar,  e_dist,   ener_res


def retarded_flux(t, direc, term_ener):
    e_mid, e_bins, MB = dist_calc(e_dist, ener_res, e_bar)
    idx = (np.abs(e_mid - (term_ener[0,0] - ener_res*0.5))).argmin()
    frac = np.sum(e_bins[idx:])
    
    """These next few lines calculate the energy flux at pellet surface"""
    flux = 0.25*dens_plas*vrms_e
    total_imp_ener = np.sum(e_bins[idx:]*e_mid[idx:])
    ener_flux = flux*total_imp_ener

    """Now to calculate the expected lifetime of such a pellet"""
    pel_num_dens = solid_dens/(m_p)
    pel_num_dens *= 10**3 # conversion from cgs to SI
    part_number = pel_num_dens*rp[t_static]**3 * 10**(-9) #conversion factor thrown in at end
    total_bond_energy = part_number*bond_energy
    ener_rate = ener_flux*(4*np.pi*rp[t_static]**2)
    ener_rate *= 10**(-6)

    lifetime = total_bond_energy/ener_rate
    return frac, ener_flux, lifetime



    