#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import glob


truth = np.loadtxt('truth.1.0')
N = 32
stoch={}
for i in xrange(N):
    filename='stoch.1.'+'{:02}'.format(i)
    print filename
    stoch[i] = np.loadtxt(filename)

assim={}
for i in xrange(N):
    filename='assim.1.'+'{:02}'.format(i)
    print filename
    assim[i] = np.loadtxt(filename)


xp = 0
for i in xrange(N):
    stoch_handle, = plt.plot(stoch[i][:,xp],color='c')
for i in xrange(N):
    assim_handle, = plt.plot(assim[i][:,xp],color='r')
for i in xrange(N):
    truth_handle, = plt.plot(truth[:,xp],color='k')

plt.legend([stoch_handle,assim_handle,truth_handle], ['Stochastic ensemble', 'Assimilation', 'Truth'],loc='best')
plt.xlabel('Timestep')
plt.ylabel('Value')
plt.show()
