#!/usr/bin/env python3

import matplotlib.pyplot as plt
from subprocess import run
import numpy as np
from matplotlib.animation import FuncAnimation

fig, _ = plt.subplots(2, 1)

points = 2 ** 16
a = 10000

X = np.linspace(0, a, points)
d12All = np.linspace(0.0008, 0.0012, 100)
D11 = np.array(points)
D12 = np.array(points)
D22 = np.array(points)
N1 = []
N2 = []
Iter = []
y11 = []
y12 = []
y21 = []
y22 = []

ax1 = plt.subplot(221)
l11, = ax1.plot(X, np.zeros(X.shape[0]), label=r'$C_{1, 1}$', animated=True)
l12, = ax1.plot(X, np.zeros(X.shape[0]), label=r'$C_{1, 2}$', animated=True)
l13, = ax1.plot(X, np.zeros(X.shape[0]), label=r'$C_{2, 2}$', animated=True)
ax1.set_xlabel(r'$d_{1, 2}$')
ax1.legend()

ax2 = plt.subplot(223)
l21, = ax2.plot([], [], label='N1', animated=True)
l22, = ax2.plot([], [], label='N2', animated=True)
ax2.set_xlabel(r'$d_{1, 2}$')
ax2.legend()

ax3 = plt.subplot(222)
l31, = ax3.plot([], [], label=r'$y_{1, 1}$', animated=True)
l32, = ax3.plot([], [], label=r'$y_{1, 2}$', animated=True)
l33, = ax3.plot([], [], label=r'$y_{2, 1}$', animated=True)
l34, = ax3.plot([], [], label=r'$y_{2, 2}$', animated=True)
ax3.set_xlabel(r'$d_{1, 2}$')
ax3.legend()

ax4 = plt.subplot(224)
l4, = ax4.plot([], [], label='iter', animated=True)
ax4.set_xlabel(r'$d_{1, 2}$')
ax4.legend()

def init():
    ax2.set_ylim([-100, 300])
    ax2.set_xlim([d12All[0], d12All[-1]])
    ax3.set_xlim([d12All[0], d12All[-1]])
    ax3.set_ylim([0.00, 0.002])
    ax4.set_xlim([d12All[0], d12All[-1]])
    ax4.set_ylim([0, 1000])
    N1 = []
    N2 = []
    Iter = []
    y11 = []
    y12 = []
    y21 = []
    y22 = []
    return l11, l12, l13, l21, l22, l31, l32, l33, l34, l4

def update(i):
    if i == 0:
        init()
    run(['./a', 'point', '-t', 'pic', '-a', str(a), '-dmat', '0.001', str(d12All[i]), '0.001', '0.001', '-sm', '0.04', '0.1', '-points', str(points)])
    D = np.genfromtxt('pic.data')
    D11 = D[0]
    D12 = D[1]
    D22 = D[2]
    l11.set_data(X, D11)
    l12.set_data(X, D12)
    l13.set_data(X, D22)
    N = np.genfromtxt('pic.N')

    N1.append(N[0])
    N2.append(N[1])

    l21.set_data(d12All[:(i+1)], N1)
    l22.set_data(d12All[:(i+1)], N2)

    O = np.genfromtxt('pic.other')
    y11.append(O[0][0])
    y12.append(O[0][1])
    y21.append(O[1][0])
    y22.append(O[1][1])
    Iter.append(O[2][0])

    l31.set_data(d12All[:(i+1)], y11)
    l32.set_data(d12All[:(i+1)], y12)
    l33.set_data(d12All[:(i+1)], y21)
    l34.set_data(d12All[:(i+1)], y22)

    l4.set_data(d12All[:(i+1)], Iter)

    y_min, y_max = ax3.get_ylim()
    if i != 0:
        ax3.set_ylim([
            min(y_min, y11[-1], y12[-1], y21[-1], y22[-1]),
            max(y_max, y11[-1], y12[-1], y21[-1], y22[-1])
            ])
    else:
        ax3.set_ylim([
            min(y11[-1], y12[-1], y21[-1], y22[-1]),
            max(y11[-1], y12[-1], y21[-1], y22[-1])
            ])
    ax1.set_ylim([
        min(np.nanmin(D11), np.nanmin(D12), np.nanmin(D22)),
        max(np.nanmax(D11), np.nanmax(D12), np.nanmax(D22))
        ])

    print('i = ' + str(i))
    return l11, l12, l13, l21, l22, l31, l32, l33, l34, l4

ani = FuncAnimation(fig, update, frames=d12All.shape[0], init_func=init, blit=True, repeat=False, interval=0)
#plt.show()
ani.save('video.mp4', dpi=400, fps=30, extra_args=['-vcodec', 'libx264'])
