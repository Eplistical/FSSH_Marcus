#!/usr/bin/env python3

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import sys


'''
mass = 2000.0
kT = 0.01
omega = 5e-4
g = 16.0
dG0 = -0.002
V = float(sys.argv[1])
'''

'''
# Landry params
mass = 1
kT = 9.5e-4
omega = 3.5e-4
dG0 = -0.015
V = np.exp(float(sys.argv[1]))
Er = 0.0239
M = np.sqrt(0.5 * Er * mass * omega**2)
g = 2 * M / mass / omega**2 
'''

# my params
mass = 1
kT = 9.5e-4
omega = 3.5e-4
dG0 = -0.015
V = np.exp(float(sys.argv[1]))
Er = 0.0239
M = np.sqrt(0.5 * Er * mass * omega**2)
g = 2 * M / mass / omega**2 


xb = 0.5 * g + dG0 / mass / omega**2 / g
Eb = 0.125 * mass * omega**2 * g**2 + 0.5 * dG0 + dG0**2 / 2 / mass / omega**2 / g**2


sigmax = np.sqrt(kT / mass / omega**2)
sigmav = np.sqrt(kT / mass)


print('V = %.6f' % (V))
print('xb = %.6f' % (xb))
print('Eb = %.6f = %.6f kT' % (Eb, Eb / kT))
print('Er = %.6f' % Er)
print('sigmax = %.6f  sigmav= %.6f' % (sigmax, sigmav))


# H00 = 0.5 * mw2 * (x + 0.5g)**2 - 0.25Er
# H11 = 0.5 * mw2 * (x - 0.5g)**2 - 0.25Er - dG0
k_marcus = V**2 * np.sqrt(np.pi / kT / Er) * np.exp(-(Er + dG0)**2 / 4 / kT / Er)

'''
print(4 * Er * kT)
print((Er + dG0)**2)
print(np.exp(-(Er + dG0)**2 / 4 / kT / Er))
print(2*np.pi / np.sqrt(4 * np.pi * Er * kT))
'''

print('log10(V) = %.6f' % (np.log10(V)))
print('ln(V) = %.6f' % (np.log(V)))
print('log10(k_marcus) = %.6f' % (np.log10(k_marcus)))
print('ln(k_marcus) = %.6f' % (np.log(k_marcus)))



def cal_H(x):
    H = np.zeros([2,2])
    H[0,0] = 0.5 * mass * omega**2 * x**2
    H[1,1] = 0.5 * mass * omega**2 * (x-g)**2 + dG0
    H[0,1] = V
    H[1,0] = np.conj(H[0,1])
    return H


Nx = 1000
xx = np.linspace(-1000,1500,Nx)
H00 = np.zeros(Nx)
H11 = np.zeros(Nx)
E0 = np.zeros(Nx)
E1 = np.zeros(Nx)

for i, x in enumerate(xx):
    H = cal_H(x)
    H00[i] = H[0,0]
    H11[i] = H[1,1]

    eva, evt = la.eigh(H)
    E0[i] = eva[0]
    E1[i] = eva[1]


plt.plot(xx, H00, '-', label='H00')
plt.plot(xx, H11, '-', label='H11')
plt.plot(xx, E0, '--', label='E0')
plt.plot(xx, E1, '--', label='E1')

plt.legend()
plt.show()
