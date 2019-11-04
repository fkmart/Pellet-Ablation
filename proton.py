#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 09:26:08 2019

@author: kyle
"""
import numpy as np
import scipy.integrate as spit


#RMEE = 0.51100334
RME = 938.2720813

e = 1.6*10**(-19)
eps0 = 8.854*10**(-12) 
#I = 21.8

M_fac = 10.0**6
KE_bot = 50
KE_top = 30000

#Distributin functino calculation

ener_res = 100
e_dist = np.arange(KE_bot, KE_top, ener_res)
e_bar = 10.0**3


def dist_calc(e_dist,ener_res, e_bar):
    MB_init = 2.07*(e_dist**(0.5)/(e_bar**1.5))*np.exp(-1.5*e_dist/e_bar)
    integ = spit.simps(MB_init, e_dist, axis = 0)
    MB_norm = MB_init/integ 
    initial_MB = np.append(e_dist, MB_norm)
    initial_MB = np.reshape(initial_MB, (int(len(initial_MB)/2),2), 'F')
    np.savetxt('electron_MB_init_norm_'+str(KE_bot)+'_'+ str(e_bar)+ '_'+str(KE_top)+'.txt', initial_MB)
    
    #Must now define the energy bins and the fraction of the distribution functino associated with that
    
    e_mid = np.zeros(len(e_dist) - 1)
    e_bins = np.zeros(len(e_dist) - 1)
    for i in range(0,len(e_mid)):
        e_mid[i] = 0.5*(e_dist[i] + e_dist[i+1])
        e_bins[i] = spit.simps((MB_norm[i], MB_norm[i+1]), (e_dist[i],e_dist[i+1]))
    return e_mid, e_bins, MB_norm
        
    #SHOULD MAKE SURE THE STOPPING CODE ACTS ON THE MID-POINT ENERGIES FOR
    #CONSISTENCY
def terminal_eedf(r, MB_norm, le_mid, particle, i):
    crit_ind = len(r) -1 
    e_term = np.zeros(le_mid)
    dist = np.flip(MB_norm.copy(), axis=0)
    e_term = []
    e_term_final = []
    dist_final = []
    for j in range(0, le_mid):
        ener_prof = np.load('EvsR_'+particle +'_t%d_E0%d.npy' % (i,j))
        ener_length = len(ener_prof)
        if ener_length < crit_ind +1:
            pass
        else:
            e_term[j] = ener_prof[le_mid -1 -j]
        if e_term[j] != 0 :
            e_term_final = np.append(e_term_final, e_term[j]*RME_e*M_fac)
            dist_final = np.append(dist_final, dist[j])
    final = np.append(e_term_final, dist_final)
    final = np.reshape(final, (len(r),2), 'F')
    np.savetxt(particle+'terminalEEDF_ebar%d_t%d.txt' % (e_bar,i), final)
    

#def stop(T,x,I,e,RMEE,RMEP):
#    dedx = -(z*e**2/(4.0*np.pi*eps0**2))*(den/(2*T))*(np.log(4.0*T/I))