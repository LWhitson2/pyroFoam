#!/usr/bin/python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
[t, Ts, Ti, Tg] = np.loadtxt("TiAvg.out", unpack=True, skiprows=1)
plt.figure(0)
plt.plot(t, Ti, 'k-')
plt.grid(b=True, which='major', color='b', linestyle='-')
plt.show()
