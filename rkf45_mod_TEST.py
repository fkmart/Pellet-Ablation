import numpy as np 
import rkf45_mod as rkf 
import matplotlib.pyplot as plt 

#test function 
#y = exp(b*y)
b = -10 
y0 = 1.0
h = 0.01 
x = np.arange(0.0,1.0 +h, h)
def f(x,y):
    return b*y 
#analytical result 
h_real = 0.0001
x_real = np.arange(0.0, 1.0 + h_real,h_real)
y_real = y0*np.exp(b*x_real)


# time loop for RKF45
nMax = 1000

xrk = np.zeros(1)
yrk = y0*np.ones(1)
hrk = np.zeros(1)
h = 0.5
for i in range(nMax):
    xloc , yloc, h , er = rkf.rkf(xrk[-1],yrk[-1],h,f,x[-1], x[0])
    xrk = np.append(xrk,xloc)
    yrk = np.append(yrk,yloc)
    if i==0:
        hrk[i] = h
    else:
        hrk = np.append(hrk,h)
    if xrk[-1]==x[-1]:
        break

fig, ax = plt.subplots()
ax.plot(x_real,y_real, color = 'teal', label = 'analytical')
ax.scatter(xrk,yrk, color = 'firebrick', label = 'rk solution')
plt.legend(loc = 'best')
plt.show()