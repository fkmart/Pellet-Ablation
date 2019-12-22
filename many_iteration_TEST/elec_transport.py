import numpy as np 
import stop_calc_rp_rc 
from gen_var import t, eps, dt
import scipy.integrate as spint 


def elec_mover(t_ind, r, eden, rden):
    rp = stop_calc_rp_rc.calc_rp(t)
    rc = stop_calc_rp_rc.calc_rc(t)
    rcd = stop_calc_rp_rc.calc_rdot(rp,rc,t)

    rp1 = rp[t_ind]
    rp2 = rp[t_ind + 1]
    rc1 = rc[t_ind]
    rc2 = rc[t_ind + 1]

    indl1 = next(p[0] for p in enumerate(r) if p[1] > rp1)
    indu1 = next(p[0] for p in enumerate(r) if p[1] > rc1) 

    indl2 = next(p[0] for p in enumerate(r) if p[1] > rp2)
    indu2 = next(p[0] for p in enumerate(r) if p[1] > rc2)
  
    dens1 = eps*((1.0 + rp1**2)/(1.0 + r[indl1:indu1]**2))
    dens2 = eps*((1.0 + rp2**2)/(1.0 + r[indl2:indu2]**2))

    dens_diff = dens1[indl1:indu1] - dens2[indl2:indu2] 

    #Now calculate the flowspeeds and fractional changes
    frac_change = np.zeros(indu1 - indl1)
    const = dens1[indu1]*rcd[t_ind]
    flowspeed = np.zeros(indu1-indl1)
    for i in range(indl1,indu1):
        frac_change[i - indl1] = dens_diff[i- indl1]/dens1[i]
        flowspeed[i - indl1] = const/dens1[i] 

    #Now interpolate the densities to all internal gridpoints 
    g = spint(rden,eden, kind = 'linear')

    ind_interp_low = next(p[0] for p in enumerate(r) if p[1] > rden[0])
    ind_interp_up = next(p[0] for p in enumerate(r) if p[1] >rden[-1])

    eden_grid = g(r[ind_interp_low:ind_interp_up])

    eden_grid_rem *= 1.0 - frac_change[:]

    eden_grid_trans = eden_grid[:]*frac_change[:]
    r_trans = r[ind_interp_low:ind_interp_up] + flowspeed[:]*dt 

    #Better way to do this is move the particles FIRST - then interpolate.
    
    

