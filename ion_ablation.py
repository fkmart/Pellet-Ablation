# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 10:28:05 2020

@author: Kyle
"""

import numpy as np 
from gen_var import diff_ion,rp_hr, N_ablat_sca, N_0_sca, delta_t, trunc_fac, I, ratio
from electron import RME, M_fac
from gen_var import v_ion, vrms_e

def ion_flux(n_ion, i, mfp_ener):
    ener_flux_n = (1.0/(RME*M_fac))*n_ion*(mfp_ener + trunc_fac*I)*(v_ion/vrms_e)
    ener_n = ener_flux_n * rp_hr[i]**2 * (1.0/ratio)
    N_A_norm = np.copy(ener_n)
    N_A_norm *= N_ablat_sca/N_0_sca
    return N_A_norm