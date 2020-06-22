import scipy.interpolate as spint 
import numpy as np
import matplotlib.pyplot as plt 

x = np.linspace(1.0, 10**5, num = 500, endpoint = 'true')

x1 = [20.0,100.0,200.0,1000.0]
D1 = [1.02, 2.41,3.63, 9.82]

f = spint.interp1d(x1,D1, kind = 'cubic', fill_value = 'extrapolate')

D = f(x)

n = 10**20 #in cm^-3
print('At T = 100K the value of D is ' + str(D[1]*10**19/n))
D_10K = D[1]*(10**5/10**2)**0.5
print('At T = 10000K the value of D is ' + str(D_10K*10**19/n))

fig,ax = plt.subplots() 

ax.plot(x1,D1, label = 'data')
ax.plot(x,D,label = 'interpolated')
plt.legend() 
plt.show()