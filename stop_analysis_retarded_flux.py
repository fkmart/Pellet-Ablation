import numpy as np
from gen_var import bond_energy, vrms_e, dens_plas, solid_dens, m_p, rp 
from gen_var import style, N_ablat_sca, N_0_sca, delta_t, tf
from electron import dist_calc, e_bar,  e_dist,   ener_res, RME, M_fac
if style == 'many':
    from gen_var import rp_hr as rp
else:
    pass

def retarded_flux(i, term_ener):
    e_mid, e_bins, MB = dist_calc(e_dist, ener_res, e_bar)
    idx = (np.abs(e_mid - (term_ener[0,0] - ener_res*0.5))).argmin()
    frac = np.sum(e_bins[idx:])
    
    """These next few lines calculate the energy flux at pellet surface"""
    #flux = 0.25*dens_plas*vrms_e
    #total_imp_ener = np.sum(e_bins[idx:]*e_mid[idx:])
    #ener_flux = flux*total_imp_ener
    ener_flux_norm = (1.0/(RME*M_fac))*np.sum(e_mid[idx:]*e_bins[idx:])
    ener_ablat_norm = ener_flux_norm*rp[i]**2
    
    """Now to calculate the expected lifetime of such a pellet"""
    #pel_num_dens = solid_dens/(m_p)
    #pel_num_dens *= 10**3 # conversion from cgs to SI
    #part_number = pel_num_dens*rp[i]**3 * 10**(-9) #conversion factor thrown in at end
    #total_bond_energy = part_number*bond_energy
    #ener_rate = ener_flux*(4*np.pi*rp[i]**2)
    #ener_rate *= 10**(-6)
    total_bonds = rp[i]**3
    ablate_rate = ener_ablat_norm*(N_ablat_sca/N_0_sca)
    lifetime = (total_bonds/ablate_rate)*delta_t*tf
    #total_bond_energy = rp[i]**3 * bond_energy/(RME*M_fac)
    #lifetime = total_bond_energy/ener_ablat_norm
    return frac, ener_ablat_norm, lifetime



    