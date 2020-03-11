import numpy as np 
import scipy.interpolate as spint 

x = np.linspace(0.0, 10.0, 100)
y = x**2 
f = spint.interp1d(x,y, kind = 'cubic')
x_new = np.asarray(2.25)
y_new = f(x_new)
print(y_new)