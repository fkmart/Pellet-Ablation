import numpy as np

def number(i):
    from gen_var import dt, dens_plas, r0, rp, tf
    from gen_var import vrms_e as vrms

    flux = 0.25*vrms*dens_plas
    if (i==0):
        no = dt*tf*flux*np.pi*r0**2
    else:
        no = dt*flux*np.pi*rp[i]**2 *r0**2
    return no