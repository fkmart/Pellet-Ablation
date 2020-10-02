# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 13:27:43 2020

@author: Kyle
"""

import numpy as np
import numba

print(numba.__version__)
print(np.__version__)
t = np.linspace(0.0,0.80,num = 100)