# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:34:34 2020

@author: Kyle
"""

import numpy as np
import os 
import csv 

direc = os.getcwd()

l = 11
r = np.linspace(0.0,10.0,num = l)
n = np.linspace(0.0,0.5, num = l)
p = np.linspace(0.5,5.5,num = l)


header = ['position','density','potential']
data = np.zeros((l,3))
data[:,0] = r
data[:,1] = n
data[:,2] = p

data_save = np.append(header,data)
data_save = np.reshape(data_save, (l+2,3))
#np.savetxt(direc + 'data_save_proper.txt', data_save, fmt = '%s', delimiter = ',')