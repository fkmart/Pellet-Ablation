import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spit

fig, ax1  = plt.subplots()
ax1.set_xlabel(r'$\tilde{t}$')
ax1.set_ylabel(r'$\tilde{r}_p$')
ax2 = ax1.twinx()
ax2.set_ylabel(r'$\tilde{r}_c$')
t = np.linspace(0.01, 0.99, 99)


def rp_calc(c):
     c1 = 1.0
     rp = c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c))))
     return rp

def rpdot(c):
     c1 = 1.0
    
     rpd = -c1*(c*0.5/(np.log(np.cosh(c))))*(np.tanh(c*(1-t)))*(1 -t)/((np.sqrt(np.log(np.cosh(c*(1-t))))))
     return rpd

def rc(c):
     r0 = np.asarray([1.0])
     c1 = 1.0
     eps = 0.01
     def rcd(t,rc):
         rcd = (eps -0.01)*((c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c)))))**2)*(-c1*(c*0.5/(np.log(np.cosh(c))))*(np.tanh(c*(1-t)))*(1 -t)/((np.sqrt(np.log(np.cosh(c*(1-t)))))))/(rc*rc)*((1.0 + 
         rc**2)/(1.0 + c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c))))**2) - 
         2.0*(c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c))))*rpd)/rc**2)*((1.0 + 
         rc**2)/(1.0 + c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c))))**2))*(rc - c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c)))) +
         np.arctan((c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c)))) - rc)/(1.0+ c1*np.sqrt((np.log(np.cosh(c*(1-t))))/(np.log(np.cosh(c))))*rc)))
         return rcd
     rc = spit.solve_ivp(rcd, (t[0], t[-1]), r0, 'RK45')
     rc = rc["y"]
     rc = rc[1:]
     return rc

c2 = [-1,-5,-10,-100]
labels = []
linestyle = ['-', '--', '-.' , ':']
count = 0
for i in c2:
     rp = rp_calc(i)
     rpd = rpdot(i)
     rc = rc(i)
     ax1.plot(t, rp, 'royalblue', linestyle[count])
     ax2.plot(t,rc, 'darkorange', linestyle[count])
     labels.append(r'$c_2 =$' + str(i))
     count +=1
plt.legend(labels)
plt.show()