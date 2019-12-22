import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as spit 

N = 9
le = 2**N + 1
print(le) 

upper = 1.5*10**4
en = np.linspace(0.0, upper, le)
e_bar = 1000.0
MB = 2.07*np.sqrt(en)*(e_bar**(-1.5))*np.exp(-1.5*en/e_bar)
print(en[1])

#define function for romberg 
def MB_func(e):
    MB = 2.07*(np.sqrt(e))*(1000.0**(-1.5))*np.exp(-1.5*e/(1000.0))
    return MB
#normalisation check 
check = spit.romb(MB, dx = en[1] - en[0])
print('Sampled integral via romb is ' + str(check))

print('evaluated value of EEDF at maximal energy is' + str(MB_func(4.0*10**4)))
check2 = spit.romberg(MB_func,0.0,upper)
print('Function integral via Romberg is ' + str(check2))
