import numpy as np
import RK4
import matplotlib.pyplot as plt


l = 500
t = np.linspace(0.0, 1.0, l)
def func(t):
    ydot = 0.0 + 0.0*t
    return ydot

y0 = 3.0

y = RK4.RK4_y(func, y0, t)

fig, ax  = plt.subplots()
ax.plot(t,y)
ax.plot(t, func(t))
plt.show()