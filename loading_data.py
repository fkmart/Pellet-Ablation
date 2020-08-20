# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:19:19 2020

@author: Kyle
"""

import os 
import numpy as np
direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration','neutral','500eV', 'analysed_outputs') + os.sep 

file = np.genfromtxt(load_dir + 'outputs5.txt', delimiter = ',',dtype = 'str')
arr1 = file[9:,1]
arr1 = np.asarray([float(i) for i in arr1])