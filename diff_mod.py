# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 10:06:47 2020

@author: Kyle
"""

import numpy as np
import scipy.signal as scisig 
from gen_var import delta_t
import matplotlib.pyplot as plt

def ion_diff(n,r,diff_ion_n,ratio):
    dt_ion = delta_t/ratio
    r_mid = r[int(len(r)*0.5)]
    Z = (r - r_mid)/(np.sqrt(dt_ion*diff_ion_n*4.0))
    f = (np.exp(-Z**2))/(np.sqrt(np.pi*diff_ion_n*dt_ion))
    n_diff = scisig.convolve(n, f, mode = 'same')
    A = np.sum(n)
    B = np.sum(n_diff)
    n_out = n_diff *A/B
    fig, ax = plt.subplots()
    ax.plot(r,n)
    ax.plot(r,f)
    ax.plot(r,n_out)

    plt.show()
    return n_diff