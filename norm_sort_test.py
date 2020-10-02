# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 08:48:10 2020

@author: Kyle
"""
import numpy as np
import matplotlib.pyplot as plt 
import scipy.interpolate as spint

fig, ax = plt.subplots()


x_new = np.linspace(0.0,8.0, num = 9, endpoint = 'true')
x_append = np.linspace(10.0,9.0, num = 2)
x = np.append(x_new, x_append)
y = x**2 - 5.0*x
ax.plot(x,y, label = 'original')

x_new = np.linspace(0.0,10.0, num = 101)

#Need to sort the initial data set from x and y 
x_sorted = np.sort(x)
y_sorted = np.zeros(len(x))
for i in range(0, len(y_sorted)):
    index = np.argmin(np.abs(x_sorted[i] - x[:]))
    y_sorted[i] = y[index]

ax.plot(x_sorted, y_sorted)
y_new = spint.pchip_interpolate( x_sorted, y_sorted, x_new)
ax.plot(x_new, y_new)

plt.show()