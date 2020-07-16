# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 16:55:08 2020

@author: Kyle
"""

def index_crit(arr):
    z = len(arr)
    for i in range(0, z):
        if arr[z - 1 - i] == 0.0:
            pass
        else:
            break
    ind_up = z - i
    for j in range(0,z):
        if arr[j] ==0.0:
            pass
        else:
            break
    if j == ind_up or j> ind_up:
        j = 0
    else: 
        pass
    ind_low = j
    return ind_up, ind_low