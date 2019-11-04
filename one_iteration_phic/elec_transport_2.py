import numpy as np 
import stop_calc_rp_rc 
from gen_var import t, eps, dt, dr
import scipy.integrate as spint 
import scipy.interpolate as spit

def elec_mover(t_ind, r, eden, rden, shift):
    rp = stop_calc_rp_rc.calc_rp(t)
    rc = stop_calc_rp_rc.calc_rc(t)
    rcd = stop_calc_rp_rc.rdot(t,rc,rp)
    lr = len(r)
   

    rp1 = rp[t_ind]
    rp2 = rp[t_ind + shift]
    rc1 = rc[t_ind]
    rc2 = rc[t_ind + shift]

    indl1 = next(p[0] for p in enumerate(r) if p[1] > rp1)
    indu1 = next(p[0] for p in enumerate(r) if p[1] > rc1) 

    indl2 = next(p[0] for p in enumerate(r) if p[1] > rp2)
    indu2 = next(p[0] for p in enumerate(r) if p[1] > rc2)
    
    dens1 = np.zeros(lr)
    dens2 = np.zeros(lr) 

    dens1[indl1:indu1] = eps*((1.0 + rp1**2)/(1.0 + r[indl1:indu1]**2))
    dens2[indl2:indu2] = eps*((1.0 + rp2**2)/(1.0 + r[indl2:indu2]**2))
    dens_diff = np.zeros(lr)
    dens_diff[indl1:indu1] = dens1[indl1:indu1] - dens2[indl1:indu1] 

    #Now calculate the flowspeeds and fractional changes
    frac_change = np.zeros(lr)
    const = dens1[indu1-1]*rcd[t_ind]
    flowspeed = np.zeros(lr)
    for i in range(indl1,indu1):
        frac_change[i] = dens_diff[i]/dens1[i]
        flowspeed[i] = const/dens1[i] 

    #Now interpolate the densities to all internal gridpoints 
    g = spit.interp1d(rden + rp1,eden, kind = 'bicubic')
    #g = spit.PchipInterpolator(rden + rp1, eden)

    ind_interp_low = next(p[0] for p in enumerate(r) if p[1] >= rden[0] + rp1)
    ind_interp_up = next(p[0] for p in enumerate(r) if p[1] >= rden[-1] + rp1) 

    eden_grid = g(r[ind_interp_low:ind_interp_up])

    ind_low = ind_interp_low - indl1 
    ind_up = indu1 -  ind_interp_up 

    eden_grid_rem = (1.0 - frac_change[ind_interp_low:ind_interp_up])*eden_grid[:]

    eden_grid_trans = eden_grid[:]*frac_change[ind_interp_low: ind_interp_up]

    r_trans = r[ind_interp_low:ind_interp_up] + flowspeed[ind_interp_low:ind_interp_up]*dt*shift

    #Better way to do this is move the particles FIRST - then interpolate.
    
    #interpolate transported values onto the same grid as for the stationary electrons
    f = spit.interp1d(r_trans, eden_grid_trans, kind = 'bicubic')
    #f = spit.pchip_interpolate(r_trans, eden_grid_trans)
    
    ind_low_2 = next(p[0] for p in enumerate(r) if p[1] >= r_trans[0])
    ind_up_2 = next(p[0] for p in enumerate(r) if p[1] >= r_trans[-1])

    eden_trans_full = f(r[ind_low_2:ind_up_2])

    #Now must define a grid from rp2 to rc2 with new charge density after transport

    #These indices are a problem!!!!!!!
    #Define in universal terms then delete zero indices
    r_new = r[indl2:indu2]
    charge_dens = np.zeros(len(r))
    charge_dens[ind_interp_low:ind_interp_up] = eden_grid_rem[:]
    charge_dens[ind_low_2:ind_up_2] += eden_trans_full[:]
    charge_dens = charge_dens[indl2:]
    charge_dens = charge_dens[:indu2-indl2]
    return charge_dens, r_new, eden_grid, r[ind_interp_low:ind_interp_up]