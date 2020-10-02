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

y2 = 3.5*x**2 + err
y2l = np.log10(y2[1:])
xl = np.log10(x)
fit = np.polyfit(x,y,deg = 1)
fit2 = np.polyfit(xl[1:],y2l,deg = 1)

print('The fitted coefficients are ' + str(fit))
print('The second fitted coefficients are ' + str(fit2))
line = fit[0]*x + fit[1]
line2 = (10**fit2[1])*x**(fit2[0]) 
fig, ax = plt.subplots()
ax.scatter(x,y)
ax.scatter(x,y2)
ax.plot(x,line)
ax.plot(x,line2)
ax.plot()
plt.show()