import charge_dens
from gen_var import rc, rp, dr
from electron import RME, M_fac
import numpy as np
import TEST, iterative_sol
"""
def elec_phi(r, q_dens, bc_in, bc_out):
    phi_0 = RME*M_fac
    
    bc_in /=phi_0
    bc_out /= phi_0
    l = len(r)
    phi = np.zeros(l)
    
    A = np.zeros((l,l))
            # setting up discretisation matrix
    for i in range(0,l):
        A[i,i] = -2.0
        #this could be tidier but not essential to change
    for i in range(1,l):
        A[i,i-1] = 1.0     
    for i in range(0,l-1):
        A[i,i+1] = 1.0  
    phi, phi_mat = TEST.new(A, phi, q_dens, bc_in, bc_out)
    return phi, phi_mat
"""
def elec_phi_test(r, q_dens):
    l = len(r)
    phi = np.zeros(l)
    phi[0] = 1.0
    phi[-1] = 1.0
    A = np.zeros((l,l))
            # setting up discretisation matrix
    for i in range(0,l):
        A[i,i] = -2.0
        #this could be tidier but not essential to change
    for i in range(1,l):
        A[i,i-1] = 1.0     
    for i in range(0,l-1):
        A[i,i+1] = 1.0  
    phi, phi_mat = iterative_sol.SOR(A, phi, q_dens, r)
    return phi, phi_mat

