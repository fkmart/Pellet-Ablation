# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 13:42:52 2020

@author: Kyle
"""
import numpy as np

def solver(n_in, D, h, delta_t):
    n_int = np.zeros(len(n_in) + 2)
    n_int[1:-1] = np.copy(n_in[:])
    D_int = np.zeros(len(n_int))
    D_int[1:-1] = D[:]
    n_new = np.zeros(len(n_int))
    for i in range(len(n_int) - 3, 0, -1):
        term1 = D_int[i] *(n_int[i+1] + n_int[i-1] - 2.0*n_int[i])/(4.0*h*h)
        term2 = ((D_int[i+1] - D_int[i-1])/(2.0*h))*(n_int[i+1] - n_int[i-1])/(2.0*h)
        n_new[i] = n_int[i] + delta_t *(term1 + term2)
    return n_new 