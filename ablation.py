
def fluxes(i,faux_dens, term_ener):
    from electron import e_dist, ener_res, e_bar, dist_calc, RME, M_fac
    import numpy as np
    from gen_var import tf, delta_t, rp_hr, N_ablat_sca, N_0_sca, e, m_e, vrms_e
    #frac = faux_dens[-1,0]
    e_mid, e_bins, MB = dist_calc(e_dist, ener_res, e_bar)
    idx = (np.abs(e_mid - (term_ener[0,0] - ener_res*0.5))).argmin()
    frac = np.sum(e_bins[idx:])
    #calculate non-dim velocities for each energy - MAY NEED REMOVED
    #v_e_nd = np.sqrt(e_mid[idx:]*e * 2.0/(m_e))
    #v_e_nd /= vrms_e
    "NORMALISED VALUES - GET REAL VALUES FROM SCALED QUANTITIES IN GEN_VAR"
    ener_flux_norm = (1.0/(RME*M_fac))*np.sum(e_mid[idx:]*e_bins[idx:])
    ener_ablat_norm = ener_flux_norm*rp_hr[i]**2
    N_A_norm = np.copy(ener_ablat_norm)
    N_A_norm *= (N_ablat_sca/N_0_sca)
    lifetime = rp_hr[i]**3 / N_A_norm
    return N_A_norm, ener_ablat_norm, lifetime, frac

def new_rp(N_A_norm,N_now_norm,i):
    from gen_var import N_0_sca, N_ablat_sca, pel_dens_numb, r0
    import numpy as np
    #N_ablat_norm = N_A_norm*(N_ablat_sca/N_0_sca)
    N_new_norm = N_now_norm - N_A_norm
    #rp_new_norm = ((3.0/(np.pi*4.0))*N_new_norm)**(1.0/3.0)
    #rp_new = (1.0/r0)*((3.0/(4.0*np.pi))*N_new_norm*N_0_sca/pel_dens_numb)**(1.0/3.0)
    rp_new = N_new_norm**(1.0/3.0)
    return rp_new, N_new_norm
    