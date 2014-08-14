#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
[t, alpha, Ts, Ti, Tg, p, EMg, Lg, Ls] = np.loadtxt("TiAvg.out", unpack=True, skiprows=1)

plt.figure(0)
plt.plot(t, Ti, 'k-', label='Ti')
#plt.plot(t, Ts, 'r-', label='Ts')
#plt.plot(t, Tg, 'b-', label='Tg')
plt.legend(loc=0)
plt.grid(b=True, which='major', color='b', linestyle='-')
plt.ylim([600, 850])
#plt.twinx()
#plt.plot(t, alpha, 'g-', label='alpha')
#plt.ylim([0, 1])

#plt.figure(1)
#plt.plot(t[:-1], (Ti[:-1] - Ti[1:])/(t[:-1] - t[1:]))
#plt.ylim([0, 10000])
#plt.grid(b=True, which='major', color='b', linestyle='-')
plt.show()
