# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:13:06 2020

@author: Kyle
"""

def writer(i,N_loss_e, N_loss_i,dens_e, dens_i,phi, stop_point, term_ener):
    from gen_var import rgl,r_grid, t_low_hr, t_upper_hr, delta_t, dt_hr, tf, rc_hr, rp_hr
    from gen_var import e_mid
    import numpy as np
    from electron import ener_res, e_bar, KE_bot, KE_top, e_dist
    e_len = len(e_mid)
    data = np.zeros((rgl,7))
    data[:,0] = r_grid
    data[:,1] = dens_e
    data[:,2] = dens_i
    data[:,3] = phi
    data[:e_len,4] = e_mid
    data[:e_len,5] = stop_point
    data[:len(term_ener),6] = term_ener
    #headers
    headers = []
    headers += ['t_low','t_up','dt', 'delta_t', 'i', 'rp','rc']
    headers += [str(t_low_hr),str(t_upper_hr),str(dt_hr), str(delta_t), str(i), str(rp_hr[i]), str(rc_hr[i])]
    headers += ['E_l', 'E_u','delta_E','E_bar','N_le', 'N_li', '']
    headers += [str(KE_bot), str(KE_top), str(ener_res), str(e_bar), str(N_loss_e), str(N_loss_i), '']
    headers += ['position', 'elec_density','ion_density','potential','initial energy','stop_point', 'stopping energy']
    headers = np.reshape(headers, (int(len(headers)/7),7))
    
    data_out = np.append(headers, data)
    data_out = np.reshape(data_out, (int(len(data_out)/7),7), 'C')
    return data_out