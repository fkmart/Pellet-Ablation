import numpy as np 
import scipy.interpolate as spint 
from gen_var import rp,rc,eps

def outer_transport(r1,r2, trans_r, k):
    rp1 = r1[0]
    rp2 = r2[0]
    rc1 = r1[-1]
    rc2 = r2[-1]

    ind_cloud = next(p[0] for p in enumerate(r2) if p[1] > rc1)

    r_range = np.linspace(r2[ind_cloud], rc2, num = len(trans_r), endpoint = 'true')
    f = spint.interp1d(trans_r, r1, kind = 'cubic')
    pre_transit_r = f(r_range)

    pre_trans_dens = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)
    d = eps*(1.0 + rp2**2)/(1.0 + r_range**2)
    d_exp = 0.01*np.exp(-k*(r_range - rp2))
    f_trans = (d - d_exp)/pre_trans_dens
    f_rem = 1.0 - f_trans 
    return f_trans, pre_transit_r

def common_transport(r1,r2,pre_transit_r,k,f_rem):
    rp1 = r1[0]
    rp2 = r2[0] 

    rc1 = r1[-1]
    rc2 = r2[-1]
    ind_int = next(p[0] for p in enumerate(r2) if p[1] >pre_transit_r[0])
    g = spint.interp1d(r2,pre_transit_r, kind = 'cubic')
    pre_transit_r2 = g(pre_transit_r)
    pre_trans_dens2 = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r2**2)
    d2 = eps*(1.0 + rp2**2)/(1.0 + pre_transit_r**2)
    d1 = eps*(1.0 + rp1**2)/(1.0 + pre_transit_r**2)
    d_exp = eps*np.exp(-k*(pre_transit_r - rp2))
    f_trans_new = (d2 - d_exp - f_rem*d1)/pre_trans_dens2 
    return f_trans_new, pre_transit_r2

def transport(r1,r2,trans_r,k, f_rem, index_low, index_up):
    rp1, rc1 = r1[0], r1[-1]
    rp2,rc2 = r2[0], r2[-1]
    r_imp = r2[index_low:index_up]
    "Need to define all the arrays that get added together in the most general equation"

    #rho_new = f_rem*rho_old + f_trans*rho_old-trans + rho_exp

    dens_new = eps*(1.0 + rp2**2)/(1.0 + r_imp**2)
    dens_old = np.zeros(np.shape(r_imp))
    dens_exp = eps*np.exp(-k*(r_imp - rp2)**2)
    pre_dens = np.zeros(np.shape(r_imp))
    if r_imp[0] > rc1 or r_imp[-1] < rp1:
        pass
    else:
        dens_old = eps*(1.0 + rp1**2)/(1.0 + r_imp**2)
    
    #need to bypass this if trans_r sits behind first point.
    g = spint.interp1d( trans_r[index_low:index_up], r1[index_low:index_up],kind = 'cubic', fill_value = 'extrapolate')
    pre_r = g(r_imp)
    pre_dens = eps*(1.0 + rp1**2)/(1.0 + pre_r**2)
    f_trans = (dens_new - dens_exp - f_rem*dens_old)/pre_dens
       
    return f_trans, pre_r