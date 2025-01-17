#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Script for analysing (basically) the structure in the flux emergence simulations
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.ndimage import gaussian_filter1d
#from fltrace import trace_fieldlines

#matplotlib.rcParams['text.usetex'] = True

#Doing (basic) diagnostics from the MF emergence runs
#Will try to automate and get parameters etc.

plot_set_number = 0  #Do this now to make things neat. Yes. Good.
runs = [0]
nsets =  len(runs)

cs = plt.cm.plasma(np.linspace(0.1,0.9,nsets))

for run in runs:

    try:
        paras = np.loadtxt('parameters/variables%03d.txt' % run)   #variables numbered based on run number (up to 1000)
    except:
        print('Run not found')
        continue

    fig_width = 15#1.0*(513.11743/72)

    nx = int(paras[1])
    ny = int(paras[2])
    nz = int(paras[3])

    nsnaps = int(paras[5])
    ndiags = int(paras[6])

    hflag = int(paras[22])

    tmax = paras[4]

    xs = np.linspace(paras[12],paras[13], nx+1)
    ys = np.linspace(paras[14],paras[15], ny+1)
    zs = np.linspace(paras[16],paras[17], nz+1)

    xc = 0.5*(xs[1:] + xs[:-1])
    yc = 0.5*(ys[1:] + ys[:-1])
    zc = 0.5*(zs[1:] + zs[:-1])

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    dz = zs[1] - zs[0]

    xc = 0.5*(xs[1:] + xs[:-1])
    yc = 0.5*(ys[1:] + ys[:-1])
    zc = 0.5*(zs[1:] + zs[:-1])

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    dz = zs[1] - zs[0]

    photo_height = 10.0
    corona_height = 20.0

    z_photo = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))
    z_corona = int((nz)*(photo_height - zs[0])/(zs[-1] - zs[0]))

    if not os.path.exists('./analysis/'):
        os.mkdir('./analysis/')

    class Grid():
        def __init__(self):
            self.x0 = xs[0]; self.x1 = xs[-1]
            self.y0 = ys[0]; self.y1 = ys[-1]
            self.z0 = zs[0]; self.z1 = zs[-1]
            self.nx = nx ; self.ny = ny; self.nz = nz

    if hflag == 0:
        data_source = '/extra/tmp/mf3d/%03d/' % run
    diag_source = './diagnostics/run%02d.nc' % run

    try:
        data = netcdf_file(diag_source, 'r', mmap=False)
    except:
        print('Failed to find data', diag_source)
        continue

    diag_titles = ['openflux', 'sumcurrent', 'avgcurrent', 'energy', 'ropeheight','avglorentz']

    ndiags = len(diag_titles)
    nrows = int(np.sqrt(ndiags-1) + 1); ncols = int((ndiags-1)/nrows) + 1
    #List diagnostics as a big array to simplify moving around the plots and things
    #Could theoretically integrate FLH stuff to this later. Maaayyybbeee.
    if plot_set_number == 0:
        #Create subplots
        fig, axs = plt.subplots(nrows, ncols)

    ts = data.variables['time'][:]
    col_num = 0; row_num = 0
    for diag in diag_titles:
        ax = axs[row_num, col_num]
        try:
            toplot = data.variables[diag][:]
        except:
            print('Data not found', diag, run)
            continue
        toplot = toplot[toplot < 1e6]

        ax.set_title(diag)
        ax.plot(ts[:len(toplot)], toplot, c = cs[plot_set_number])

        col_num = (col_num + 1)%ncols
        if col_num == 0:
            row_num += 1

        ax.set_xlim(0, tmax)

    data.close()

    if plot_set_number == nsets-1:
        plt.tight_layout()
        plt.savefig('diagnostics.png')
        plt.show()
        plt.close()

    plot_set_number += 1

if True:
    #Plot Lorentz force for given set of snapshots, jsut based on the average over a horizontal slice (no need to be fancy)
    if not os.path.exists('lorentzplots'):
        os.mkdir('lorentzplots')

    nsnaps = 399
    for snap in range(nsnaps):
        print('Snap number ', snap)
        allmax = 0.0
        for ri, run in enumerate(runs):
            if ri == 0:
                fig, axs = plt.subplots(1, 1)

            diag_source = './diagnostics/run%02d.nc' % run

            try:
                data = netcdf_file(diag_source, 'r', mmap=False)
            except:
                print('Failed to find data', diag_source)
                continue

            lf = data.variables['lfheights'][:].T

            #print(run, np.sum(lf))
            end = len(ts[ts < 1e6])
            lf = lf[:end,:]

            allmax = max(allmax, np.max(lf[:end,:]))

            plt.plot(lf[snap,:], c = cs[run], label = run)

            if ri == len(runs) - 1:
                plt.title('t = %02d' % data.variables['time'][snap])
                plt.legend()
                plt.ylim(0.0, allmax*1.2)

                plt.tight_layout()
                plt.savefig('lorentzplots/%03d.png' % snap)
                plt.close()

        if snap == nsnaps-1:
            #Do the video
            os.system('ffmpeg -y -framerate 20 -i ./lorentzplots/%03d.png -b:v 10M lorentz.mp4')

            if True:
                os.system('rm -r lorentzplots')

if False:
    #Plot Lorentz Forces (colormesh maybe?!?, Alice is sometimes wise...)
    plot_set_number = 0  #Do this now to make things neat. Yes. Good.
    row_num = 0; col_num = 0
    nplots = len(runs)
    nrows = int(np.sqrt(nplots) + 1); ncols = int((nplots-1)/nrows) + 1

    for run in runs:

        #List diagnostics as a big array to simplify moving around the plots and things
        #Could theoretically integrate FLH stuff to this later. Maaayyybbeee.

        if plot_set_number == 0:
            #Create subplots
            fig, axs = plt.subplots(nrows, ncols)

        if np.ndim(axs) > 1:
            ax = axs[row_num, col_num]
        else:
            ax = axs[row_num]

        print(nrows, ncols)
        if hflag == 0:
            data_source = '/extra/tmp/mf3d/%03d/' % run
        diag_source = './diagnostics/run%02d.nc' % run

        try:
            data = netcdf_file(diag_source, 'r', mmap=False)
        except:
            print('Failed to find data', diag_source)
            continue

        lf = data.variables['lfheights'][:].T

        #print(run, np.sum(lf))
        end = len(ts[ts < 1e6])
        lf = lf[:end,:]

        #Various variations that it might be worthwhile plotting.
        #Could normalise based on the mean at each time?

        #for n in range(end):
        #    lf[n,:] = lf[n,:]/np.max(lf[n,:])

        im = ax.pcolormesh(ts[:end], zc, lf.T)

        ax.set_title('Run number = %d' % run)
        plt.colorbar(im, ax = ax)
        data.close()

        if plot_set_number == nsets-1:
            plt.tight_layout()
            plt.show()

        col_num = (col_num + 1)%ncols
        if col_num == 0:
            row_num += 1

        plot_set_number += 1














