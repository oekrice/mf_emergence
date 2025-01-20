#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27


Script for plotting the input magnetograms -- just the ones from LARE. For  awrite-up of sorts...
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
matplotlib.rcParams['text.usetex'] = True

mag_nums = [150,300,450]

nmags = len(mag_nums)

fig, axs = plt.subplots(nmags, 3, figsize = (10,7))

x0 = -130; x1 = 130
y0 = -130; y1 = 130

for mi, mag_num in enumerate(mag_nums):

    time = mag_num*0.5

    source = '../magnetograms/%04d.nc' % mag_num

    data = netcdf_file(source, 'r', mmap=False)

    bx = data.variables['bx'][:]
    by = data.variables['by'][:]
    bz = data.variables['bz'][:]

    data.close()

    nx = bz.shape[0]; ny = bz.shape[1]
    xs = np.linspace(x0, x1, nx+1)
    ys = np.linspace(y0, y1, ny+1)

    xc = np.linspace(x0, x1, nx+2)
    yc = np.linspace(y0, y1, ny+2)

    toplot = bx
    ax = axs[mi,0]
    im = ax.pcolormesh(xc, ys, toplot.T, rasterized = True)
    plt.colorbar(im, ax=ax)
    ax.set_title('$B_x$, t = %d' % time)

    toplot = by
    ax = axs[mi,1]
    im = ax.pcolormesh(xs, yc, toplot.T, rasterized = True)
    plt.colorbar(im, ax=ax)
    ax.set_title('$B_y$, t = %d' % time)

    toplot = bz
    ax = axs[mi,2]
    im = ax.pcolormesh(xs, ys, toplot.T, rasterized = True)
    plt.colorbar(im, ax=ax)
    ax.set_title('$B_z$, t = %d' % time)


plt.tight_layout()
plt.savefig('magplots.pdf')
plt.savefig('magplots.png')
plt.show()
plt.close()

