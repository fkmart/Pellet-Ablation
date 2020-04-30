import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as spint 

x1 = np.linspace(0.0,10.0, num = 10, endpoint = 'true')
y1 = 4.0*x1 
x_new = np.asarray(8.45)

f = spint.interp1d(x1,y1,kind = 'cubic', fill_value = 'extrapolate')

y_new = f(x_new)

fig, ax  = plt.subplots() 
ax.plot(x1,y1, linestyle = '--')
ax.scatter(x_new,y_new, color = 'red')
plt.show()