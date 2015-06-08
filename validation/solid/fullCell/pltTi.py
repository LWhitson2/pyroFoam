#!/usr/bin/python
import numpy as np
import scipy as sp
import scipy.special
import matplotlib as mpl
import matplotlib.pyplot as plt
[t, alpha, Ts, Ti, Tg, p, EMg, Lg, Ls] = np.loadtxt("TiAvg.out", unpack=True, skiprows=1)

T0 = 500.
T0i = 700.
q = 1311287.64368
k = 0.2
cp = 1400.
rho = 1800.
alpha = k/(rho*cp)

loc = 0.00001

def Texact(t,x):
  # T = 2*np.sqrt(alpha*t/np.pi)*np.exp(-x**2/(4*alpha*t))
  # T += x*(sp.special.erf(x/(2*np.sqrt(alpha*t))) - 1)
  # T *= q/k
  T = sp.special.erfc(x/(2*np.sqrt(alpha*t)))
  T *= (T0i - T0)
  return T + T0

plt.figure(0)
plt.plot(t, Ti, 'k-', label='Ti')
plt.plot(t, Ts, 'r-', label='Ts')
plt.plot(t, Texact(t,0), 'b-', label='TiExact')
plt.plot(t, Texact(t,loc), 'g-', label='TsExact')
plt.legend(loc=0)
plt.grid(b=True, which='major', color='b', linestyle='-')
# plt.ylim([300, 850])
#plt.twinx()
#plt.plot(t, alpha, 'g-', label='alpha')
#plt.ylim([0, 1])

plt.figure(1)
plt.plot(t, np.abs(Ti - Texact(t,0))/(Texact(t,0) - T0), 'k-', label='Ti Error')
plt.plot(t, np.abs(Ts - Texact(t,loc))/(Texact(t,loc) - T0), 'r-', label='Ts Error')
plt.yscale('log')
plt.ylim([1.e-4, 1.])
plt.legend(loc=0)
plt.show()

# t = np.linspace(0,1,1000)
# print Texact(t, 0.00199)
# print Texact(t, 0.)