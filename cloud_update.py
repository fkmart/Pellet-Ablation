import numpy as np 

import discret_mat
import iterative_sol as SOR 
import elec_transport_2 as e_trans 
import common_interpolator as com_int 

elec_interp, ind_low, ind_up = com_int.common_interp(r_internal, real_dens[:,1], real_dens[:,0]) # interpolating charge onto all gridpoints
acc_elec_dens[ind_low:ind_up] = elec_interp[:] + acc_elec_dens[ind_low:ind_up]

import cloud_mover 

low = next(p[0] for p in enumerate(r) if p[1] > rp[j])
up = next(p[0] for p in enumerate(r) if p[1] > rc[j])

A = discret_mat.discret(r[ind1_new:ind2_new]) #r here needs to be from just in front of the pellet to the cloud at the new time.
phic = SOR.SOR(A, pot, acc_elec_dens, r[ind1_new:ind2_new])