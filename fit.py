#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib as mpl
HAVE_DISPLAY = 0 == os.system('python3 -c "import matplotlib.pyplot as plt;plt.figure()"')
if not HAVE_DISPLAY:
    mpl.use('Agg')
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, k, C):
    return a * np.exp(k*x) + C
p0 = (1,-1e-6,1e-3)

fig = plt.figure()


fname = sys.argv[1]
a = np.loadtxt(fname)
Nr = -1
tarr = a[:Nr,0]
data = a[:Nr,1]

plt.plot(tarr, data)
popt, pcov = curve_fit(func, tarr, data, p0=p0)
k = popt[1]

print('log10(k) = %16.4f' % np.log10(-k))
print('ln(k) = %16.4f' % np.log(-k))

xx = np.linspace(tarr[0], tarr[-1] , 2000)
yy = func(xx, *popt)
plt.plot(xx, yy, '--')

plt.show()
