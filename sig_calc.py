import numpy as np 
def sig_calc(A,r0,rc,A_low):
    k = (rc - r0)/(np.sqrt(2.0*np.log(A/A_low)))
    return k 
    