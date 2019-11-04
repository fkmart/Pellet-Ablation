import numpy as np 
import scipy.integrate as spinteg
import scipy.interpolate as spint 
import os

def phi_calc(t, particle, arr3):
    import SOR_lapl_1D
    import MB_calc 
    import number 
    from gen_var import rc,rp, phi_plas, phi_p, dr, epsilon0, r0
    from electron import RME, M_fac, e_dist, e_bar

    """Setting up essential variables"""
    rcloud = rc[t]
    rpellet = rp[t]
    r = np.arange(0,rcloud-rpellet,dr)
    l = len(r)
    rel_tol = 1e-6
   

    """Takes the stopping points and their associated energies, interpolates the density profile over the finer grid"""
    stop_points = arr3[0:]
    s_p = np.zeros(len(stop_points))
    for k in range(0, len(s_p)):
        s_p[k] = stop_points[k][0]
    stop_points = s_p
    mb_dist = MB_calc.MB(e_dist,e_bar)
    mb_dist_flip = np.flip(mb_dist, axis = 0)
    den_cloud =  spint.pchip_interpolate(np.flip(stop_points, axis = 0),mb_dist_flip[-len(stop_points):], r)

    """Determining flux of particles into the domain"""

    no = number.number(t)

    """Quick check to determine the fractin fo particles from the background that stops in the cloud"""

    stop_frac = spinteg.simps(den_cloud[:len(stop_points)],e_dist[:len(stop_points)])

    """Charge density calculated and normalised"""

    q_dens = no*den_cloud
    tot = spinteg.simps(q_dens, r)
    q_dens /= tot

    """SOR method applied"""
    #Getting normalised variables/setting up arrays
    phi = np.zeros(l-2)
    phi_0 = M_fac*RME 
    phi_p /= phi_0
    phi_plas /= phi_0
    rd = r0*(rcloud - rpellet)
    norm = rd**2/(epsilon0*phi_0)
    q_dens *= norm 

    phi = SOR_lapl_1D.SOR(phi, q_dens, rel_tol, phi_p, phi_plas)
    np.savetxt(os.path.join('Analysed Outputs', 'GNU_'+particle+'cloud_potential_1keV_time_'+str(t)+'.txt'), phi*phi_0)
    return phi
