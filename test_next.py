import numpy as np 

y = np.linspace(0, 10, 11)
print(y)

ob = 3.2

s = next(x[0] for x in enumerate(y) if x[1] > ob) 
print(s)