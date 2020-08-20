# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:34:34 2020

@author: Kyle
"""

import numpy as np
import os  

direc = os.getcwd()

l = 11
r = np.linspace(0.0,10.0,num = l)
n = np.linspace(0.0,0.5, num = l)
p = np.linspace(0.5,5.5,num = l)

header1 = ['t_low','t_hi','dt']
header2 = ['1e-4','1e-2','1e-5']
header3 = ['delta t','tf','N_loss']
header4 = ['1e-4','1e-3','1e19']
header5 = ['delta_E','e_bar','r_p']
header6 = ['100','1e3','0.9']
header7 = ['r_c','','']
header8 = ['9.0','','']

headers = header1 + header2 + header3 + header4 + header5 + header6 + header7 + header8
last_header = headers + ['position','density','potential']

final_header = np.reshape(last_header,(int(len(last_header)/3),3))

data = np.zeros((l,3))
data[:,0] = r
data[:,1] = n
data[:,2] = p

data_out = np.append(final_header,data)
data_out = np.reshape(data_out, (int(len(data_out)/3),3),'C')
#np.savetxt(direc + 'data_save_proper.txt', data_save, fmt = '%s', delimiter = ',')