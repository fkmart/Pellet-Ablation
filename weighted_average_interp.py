# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:58:45 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0, 1.0, num  = 101)
x_wide = np.linspace(0.0,1.0, num = 6)

y = 0.01*x
y_wide = 0.01*x_wide

#forced interpolation

x_new = 0.427

def makeshift(x_new,x,y):    
    i= np.argmin(np.abs(x - x_new))
    
    if x[i] < x_new:
        x_low = x[i]
        x_up = x[i+1]
        y_low = y[i]
        y_up = y[i+1]
    else:
        x_up = x[i]
        x_low = x[i-1]
        y_up = y[i]
        y_low = y[i-1]
    
    x_diff = np.abs(x_low - x_up)
    
    low_scale = 1.0 - np.abs((x_new - x_low)/x_diff)
    up_scale = 1.0 - np.abs((x_new - x_up)/x_diff)
    
    
    y_new = low_scale*y_low + up_scale*y_up
    return y_new

y_new = makeshift(x_new, x_wide, y_wide)

fig, ax = plt.subplots()
ax.plot(x,y)
ax.scatter(x_wide,y_wide, marker = 'x', color = 'orange')
ax.scatter(x_new, y_new, marker = 'x', color = 'violet')
plt.show()