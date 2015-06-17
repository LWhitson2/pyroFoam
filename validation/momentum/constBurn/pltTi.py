#!/usr/bin/python
import numpy as np
import scipy as sp
import scipy.special
import matplotlib as mpl
import matplotlib.pyplot as plt
[t, alpha, Ts, Ti, Tg, p, EMg, rhog, Uy] = np.loadtxt("TiAvg.out", unpack=True, skiprows=1)

T0 = 500.
T0i = 700.
q = 1311287.64368
k = 0.2
cp = 1400.
rho = 1800.
a = k/(rho*cp)
alphaMin = 0.05

loc = 0.00001

def Texact(t,x):
  # T = 2*np.sqrt(alpha*t/np.pi)*np.exp(-x**2/(4*alpha*t))
  # T += x*(sp.special.erf(x/(2*np.sqrt(alpha*t))) - 1)
  # T *= q/k
  T = sp.special.erfc(x/(2*np.sqrt(a*t)))
  T *= (T0i - T0)
  return T + T0

plt.figure(0)
plt.plot(t, p, 'k-', label='Pressure')
plt.legend(loc=0)
plt.grid(b=True, which='major', color='b', linestyle='-')
# plt.ylim([300, 850])
plt.twinx()
plt.plot(t, alpha, 'g-', label='alpha')
plt.plot(t, np.ones(len(t))*alphaMin, 'g--')
plt.ylim([0, 1])

plt.figure(1)
plt.plot(t, rhog, 'k-', label='Density')
plt.legend(loc=0)
plt.grid(b=True, which='major', color='b', linestyle='-')
# plt.ylim([300, 850])
plt.twinx()
plt.plot(t, alpha, 'g-', label='alpha')
plt.plot(t, np.ones(len(t))*alphaMin, 'g--')
plt.ylim([0, 1])

plt.figure(2)
plt.plot(t, Uy, 'k-', label='Y-Velocity')
plt.legend(loc=0)
plt.grid(b=True, which='major', color='b', linestyle='-')
# plt.ylim([300, 850])
plt.twinx()
plt.plot(t, alpha, 'g-', label='alpha')
plt.plot(t, np.ones(len(t))*alphaMin, 'g--')
plt.ylim([0, 1])
plt.show()