import numpy as np 
import scipy.integrate as spint 
import matplotlib.pyplot as plt

lx = 65
x = np.linspace(0.0, 32, lx)
y = np.zeros(lx)
y[0:33] = x[0:33]**2

print('Analytically, the integral of x**2 between 0.0, and 16.0 is '+ str((16.0**3)/3.0 ))
I_simps = spint.simps(y,x,dx = x[1] - x[0])
I_romb = spint.romb(y,dx = x[1] - x[0])
print('Numerical integration by Simpson yields ' + str(I_simps))
print('Numerical integration by Romberg yields ' + str(I_romb))

fig, ax  = plt.subplots()
ax.plot(x,y)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$', rotation = 0)
plt.show()