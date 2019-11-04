import numpy as np
import scipy.integrate as spint

def MB(E, e_bar):
    f = (2.0/(np.sqrt(np.pi)))*(3.0/(2.0*e_bar))**(1.5)*E**0.5 * np.exp(-((3.0*E)/(2.0*e_bar)))
    I = spint.simps(f,E)
    f /= I
    return f
