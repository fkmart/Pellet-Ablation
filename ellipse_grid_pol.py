import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.tri as mtri

#ellipse params 

a = 1.0
b = 2.0

rad = np.linspace(-b, b, endpoint = 'true', num = 51)
theta = np.linspace(0.0, 2*np.pi, endpoint = 'true', num = 51)

R,T = np.meshgrid(rad,theta)
#X,Y = X.flatten(), Y.flatten()
ell = (b*a)/np.sqrt((b*np.cos(T))**2 + (a*np.sin(T))**2)

"""
#Define the colourmesh 
tri = mtri.Triangulation(X, Y)
fig = plt.figure()
plt.figaspect(a/b)
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_trisurf(X, Y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
el = b*np.sqrt(1.0 - (x/a)**2)
ax.plot(x,el,np.ones(len(x)), color = 'black', linestyle = '--')
plt.show()"""

fig, ax  = plt.subplots(subplot_kw={'projection': 'polar'}) 
ax.pcolormesh(T,R,ell, cmap = 'viridis')
ax.set_yticks([])
ax.set_xticks([])
plt.show()

