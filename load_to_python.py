# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:38:45 2020

@author: Kyle
"""

import numpy as np 
import os 

direc = os.getcwd()
load_dir = os.path.join(direc,'many_iteration','neutral', '1000eV','analysed_outputs') + os.sep
file1 = np.loadtxt(load_dir + 't_params.txt', delimiter = ',', dtype = 'str')

print(file1[0])
file = np.loadtxt('savetest.txt', dtype = 'str', delimiter = ',')
number = file[1,:]
numbers = [float(i) for i in number]
