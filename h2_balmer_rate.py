# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 11:49:59 2020

@author: Kyle
"""
import numpy as np
import scipy.interpolate as spint

def balmer_rate(e_eval):
    file_eff = np.loadtxt('h2_effective.txt')
    e_eff = file_eff[:,0]
    cs_eff = file_eff[:,1]
    cs_eval_eff = spint.pchip_interpolate(e_eff,cs_eff,e_eval)
    
    file_ion = np.loadtxt('h2_ionisation.txt')
    e_ion = file_ion[:,0]
    cs_ion = file_ion[:,1]
    cs_eval_ion = spint.pchip_interpolate(e_ion,cs_ion, e_eval)
    
    file_elas = np.loadtxt('h2_elastic.txt')
    e_elas = file_elas[:,0]
    cs_elas = file_elas[:,1]
    cs_eval_elas = spint.pchip_interpolate(e_elas, cs_elas, e_eval)
    
    cs_eval_balm = 0.0
    i = 15
    file = np.loadtxt('h2_excitation' + str(i) + '.txt')
    ener = file[:,0]
    cs = file[:,1]
    cs_eval_balm += spint.pchip_interpolate(ener,cs, e_eval) 
    return cs_eval_balm/cs_eval_eff