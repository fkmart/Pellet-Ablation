# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:34:53 2020

@author: Kyle
"""

import numpy as np

text = ['time (s)','space (m)','speed (m/s)']
numbers = ['5.0','50.0','10.0']

savefile = [text,numbers]
np.savetxt('savetest.txt',savefile, delimiter = ',', fmt = '%s')