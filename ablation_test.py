# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:08:34 2020

@author: Kyle
"""

import os
import numpy as np 
import ablation as ab
from gen_var import rp, N_0_sca, pel_dens_numb, r0

direc = os.getcwd() 
load_dir = os.path.join(direc, 'many_iteration', 'poisson','analysed_outputs','t_20') + os.sep
faux_dens = np.loadtxt(load_dir + 'density_t20.txt')
term_ener = np.loadtxt(load_dir + 'terminal_energy_t20.txt')


N_now = (4.0/3.0)*np.pi*rp[20]**3 * pel_dens_numb * r0**3
N_now_n = N_now/N_0_sca
N_An = ab.fluxes(20,faux_dens, term_ener)
rpn = ab.new_rp(N_An, N_now_n, 20)