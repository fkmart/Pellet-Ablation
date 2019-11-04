import numpy as np

from gen_var import *
from electron import RME as RME_e
from proton import *

#Use a very simple test case to guarantee that you get the correct number out in the end. Calculate explicitly and normalise agian afterwards if need be.

e = 1.6*10**(-19)
mp = 1.67*10**(-27)
me = 9.11*10**(-31)


E = 20000
E_n = E/(RME*10**6) 
E_r = E*e 

v = np.sqrt(2.0*E_r/mp)

logarg = (me*v**2)/(I*e)
x = 'done'

#Now try to get the same answer in normalised units

numer = 2.0*E_n*RME_e/RME
denom = I/(RME*M_fac)
logarg2 = numer/denom
x = 8
