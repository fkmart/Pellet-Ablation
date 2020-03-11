import numpy as np 
from gen_var import rc,rp , eps, dt
import stop_calc_rp_rc
from gen_var import t as t_arr

def e_mover(t, shift, grid1,den_elec):
    rcd, other = stop_calc_rp_rc.rcdot(t_arr,rc,rp)

    rp1 = rp[t] #cloud and pellet sizes at each time
    rp2 = rp[t + shift]
    rc1 = rc[t]
    rc2 = rc[t + shift]

    grid2 = np.linspace(rp2,rc2, num = len(grid1), endpoint = 'true')

    n_dens1 = eps*((1.0 + rp1**2)/(1.0 + grid1[:]**2))

    const = n_dens1[-1]*rcd[0,t]
    rdot = np.zeros(len(grid1))
    frac_change = (1.0 + rp2**2)/(1.0 + rp1**2)
    for i in range(0, len(grid1)):
        rdot[i] = const/n_dens1[i]

    elec_rem = frac_change*den_elec
    elec_tran = den_elec - elec_rem 

    r_trans = grid1 + rdot*dt*shift 
    #need to append and reshape remaining and tranpsorted arrays
    rem_arr = np.append(elec_rem, grid1)
    rem_arr = np.reshape(rem_arr, (int(0.5*len(rem_arr)),2), order = 'F')
    tran_arr = np.append(elec_tran, r_trans)
    tran_arr = np.reshape(tran_arr, (int(0.5*len(tran_arr)),2), order = 'F')
    return rem_arr, tran_arr
