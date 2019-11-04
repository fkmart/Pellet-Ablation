import numpy as np 

def cloud_mover(e_dens,i,shift):
    from gen_var import rp, rc 

    rp1 = rp[i] 
    rp2 = rp[i+shift]
    rc1 = rc[i] 
    rc2 = rc[i+shift]