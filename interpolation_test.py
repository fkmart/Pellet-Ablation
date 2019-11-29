import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spit
import os 
from gen_var import t_start, rp, rc, pel_pot, lp, sig

x = np.linspace(0.0,1.0, 10)
y =  -x**3 + x**2 + 3.0

f = spit.interp1d(x,y,kind = 'linear')
x_new = np.linspace(0.0, 1.0, 100)
y_new = f(x_new)

g = spit.interp1d(x,y, kind = 'cubic')
y2_new = g(x_new)

pchip_func = spit.PchipInterpolator(x,y, extrapolate = 'true')
y_pchip = pchip_func(x_new)

fig, ax = plt.subplots()
ax.scatter(x,y, label = 'original', marker = 'x')
ax.plot(x_new, y_new, label = 'linear')
#ax.plot(x_new,y2_new, label = 'cubic' )
ax.plot(x_new, y_pchip, label = 'p_chip')
plt.legend()
plt.show()