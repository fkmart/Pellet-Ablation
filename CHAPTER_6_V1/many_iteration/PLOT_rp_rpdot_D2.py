import numpy as np 
import matplotlib.pyplot as plt 
import math as mt 

k2 = -10.0 
k = 1.0 

t = np.linspace(0.0, 1.0, 1000)

rp_d2 = k*np.sqrt(1 - t)
rpd_d2 = -0.5*k*1.0/(np.sqrt(1.0 - t))

rp = 