import numpy as np
from gen_var import rc
s = np.linspace(0.0,3200, num = 1000)

r1 = np.linspace(0.0,rc[1], num = len(s), endpoint = 'true')
print('Cloud size is ' + str(r1[-1]))
print(s[1])
print ( r1[1])