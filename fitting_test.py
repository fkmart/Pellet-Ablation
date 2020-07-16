# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:20:05 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
import random

x = np.linspace(0,10,11)
err = np.zeros(11)
for i in range(0,11):
    err[i] = 0.1*random.randrange(1,3)
y = 2.0*x + 3.0 + err

fit = np.polyfit(x,y,deg = 1)

line = fit[0]*x + fit[1]
fig, ax = plt.subplots()
ax.scatter(x,y)
ax.plot(x,line)
plt.show()