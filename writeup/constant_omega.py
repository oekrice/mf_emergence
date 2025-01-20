#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Script for plotting the helicity analysis for constant omega runs. Data assumed to be copied from Hamilton, which can be large.

Don't need to do rope height etc.
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file

from scipy.interpolate import RegularGridInterpolator
from scipy.fft import fft, ifft2, fft2, ifft

matplotlib.rcParams['text.usetex'] = True

#There are seven runs. Use the same script here and on Hamilton for neatness purposes.

class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self, run):

        paras = np.loadtxt('../parameters/variables%03d.txt' % run)

        import_resolution = 128

        #Define the grid onto which the electric field should be outputted
        self.nx = import_resolution
        self.ny = import_resolution

        self.x0 = paras[12]; self.x1 = paras[13]
        self.y0 = paras[14]; self.y1 = paras[15]

        self.xs = np.linspace(self.x0, self.x1, self.nx+1)
        self.ys = np.linspace(self.y0, self.y1, self.ny+1)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]

        self.xc = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, self.nx+2)
        self.yc = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, self.ny+2)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]

        self.xs_import = np.linspace(self.x0, self.x1, import_resolution+1)
        self.ys_import = np.linspace(self.y0, self.y1, import_resolution+1)


        self.dx_import = self.xs[1] - self.xs[0]
        self.dy_import = self.ys[1] - self.ys[0]

        self.xc_import = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, import_resolution+2)
        self.yc_import = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, import_resolution+2)


class compute_inplane_helicity():
    #Uses all three components of the magnetic field to give values for A.B at EACH of the input timesteps.
    #Requires the script from Chris' flt code
    def __init__(self, run_min, run_max):

        grid = Grid(0)
        h_ref = []
        h_all = [[] for _ in range(run_max - run_min)]
        ts = []
        #Call run = -1 the reference case
        runs = [-1] + np.arange(run_min, run_max).tolist()
        omegas = []
        for ri, run in enumerate(runs):

            if run >= 0:
                paras = np.loadtxt('../parameters/variables%03d.txt' % run)
                omega = paras[28]
                omegas.append(omega)
                print(run, 'omega', omega)

            for snap in range(0,500,10):

                if (snap%10) == 0:
                    print('Run', run, ', snap', snap)
                #Find from LARE

                if run < 0:
                    source = '../magnetograms/'
                else:
                    source = '../mf_mags/%03d/' % run

                data_directory = source

                bfield_fname = '%s%04d.nc' % (data_directory, snap)

                try:
                    data = netcdf_file(bfield_fname, 'r', mmap=False)

                except:
                    print('File', bfield_fname, 'not found')
                    continue

                if not run < 0:
                    bx_in = np.swapaxes(data.variables['bx'][:],0,1)
                    by_in = np.swapaxes(data.variables['by'][:],0,1)
                    bz_in = np.swapaxes(data.variables['bz'][:],0,1)

                else:   #These ones are already swapped for some reason I decided in the past
                    bx_in = data.variables['bx'][:]
                    by_in = data.variables['by'][:]
                    bz_in = data.variables['bz'][:]
                #Trim the edges out as these can go a bit screwy and bugger up the results
                bx = np.zeros((np.shape(bx_in)))
                by = np.zeros((np.shape(by_in)))
                bz = np.zeros((np.shape(bz_in)))

                bx[1:-1,1:-1] = bx_in[1:-1,1:-1]
                by[1:-1,1:-1] = by_in[1:-1,1:-1]
                bz[1:-1,1:-1] = bz_in[1:-1,1:-1]

                #Need to average these to grid centres to get the FFT to work
                bx0 = 0.5*(bx[1:,:] + bx[1:,:])
                by0 = 0.5*(by[:,1:] + by[:,:-1])
                bz0 = bz[:,:]

                mag_nx = bz.shape[0]; mag_ny = bz.shape[1]
                mag_dx = (grid.xs[-1] - grid.xs[0])/mag_nx
                mag_dy = (grid.ys[-1] - grid.ys[0])/mag_ny


                def norm2d(vec):
                    mag = np.linalg.norm(vec)
                    if (mag > 0.0):
                        v = vec/mag
                    else:
                        v = np.array([0, 0])
                    return np.array([v[0],v[1],0.0])

                def getFrequencyMatrix(ncells,spacing):
                    freqlist1da =np.roll(np.linspace(-ncells[0]/2,ncells[0]/2-1,ncells[0]),round(ncells[0]/2))/(ncells[0]*spacing[0])
                    freqlist1db =np.roll(np.linspace(-ncells[1]/2,ncells[1]/2-1,ncells[1]),round(ncells[1]/2))/(ncells[1]*spacing[1])
                    return np.array([np.array([np.array([2.0*np.pi*freqlist1da[i],2.0*np.pi*freqlist1db[j]]) for j in range(len(freqlist1db))]) for i  in range(len(freqlist1da))]);

                #Find in -plane vector potential in the winding gauge

                fm = getFrequencyMatrix([mag_nx, mag_ny],[mag_dx, mag_dy]);
                # make the basis

                kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
                kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])
                # note in the k matrix below the k=0 element is set to one so we can divide by it.
                k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T

                nx = bz.shape[0]; ny = bz.shape[1]
                aftx = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)
                afty = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)
                aftz = np.zeros([bz.shape[0],bz.shape[1]],dtype=np.complex128)

                fbx = fft2(bx0[:,:]); fby = fft2(by0[:,:]); fbz = fft2(bz0[:,:])

                akperp = -1j*fbz/k
                ## fix i =j  element
                akw = 1j*(-(kparr[:,:,1])*fbx + (kparr[:,:,0])*fby)/k
                ## fix i =j  element
                aftx[:,:] = akperp*kperp[:,:,0]
                afty[:,:] = akperp*kperp[:,:,1]
                aftz[:,:] = akperp*kperp[:,:,2]+akw

                ax0 = ifft2(aftx[:,:])
                ay0 = ifft2(afty[:,:])
                az0 = ifft2(aftz[:,:])
                ax0 = np.real(ax0)
                ay0 = np.real(ay0)
                az0 = np.real(az0)

                ax = np.zeros((nx, ny+1))
                ay = np.zeros((nx+1, ny))

                ax[:,1:-1] = 0.5*(ax0[:,1:] + ax0[:,:-1])
                ay[1:-1,:] = 0.5*(ay0[1:,:] + ay0[:-1,:])

                bz_test = (ay[1:,:] - ay[:-1,:])/mag_dx - (ax[:,1:] - ax[:,:-1])/mag_dy

                #print('Vector potential test', np.max(np.abs(bz[1:-1,1:-1] - bz_test[1:-1,1:-1]))/np.max(np.abs(bz[1:-1,1:-1])))
                #This vector potential should be reasonably OK... Need code to test though

                hfield = ax0*bx0 + ay0*by0 + az0*bz0

                if snap == 300 and ri < 4:

                    if run < 0:
                        fig, axs = plt.subplots(4,4, figsize = (10,10))

                    if run < 0:
                        row = 0

                    else:
                        row = ri - 1

                    ax = axs[ri, 0]
                    im = ax.pcolormesh(bx.T, cmap = 'seismic', vmax = np.max(np.abs(bx)), vmin = -np.max(np.abs(bx)))
                    plt.colorbar(im, ax = ax)
                    ax = axs[ri, 1]
                    im =axs[ri,1].pcolormesh(by.T, cmap = 'seismic', vmax = np.max(np.abs(by)), vmin = -np.max(np.abs(by)))
                    plt.colorbar(im, ax = ax)
                    ax = axs[ri, 2]
                    im = axs[ri,2].pcolormesh(bz.T, cmap = 'seismic', vmax = np.max(np.abs(bz)), vmin = -np.max(np.abs(bz)))
                    plt.colorbar(im, ax = ax)
                    ax = axs[ri, 3]
                    im = axs[ri,3].pcolormesh(hfield.T, cmap = 'seismic', vmax = np.max(np.abs(hfield)), vmin = -np.max(np.abs(hfield)))
                    plt.colorbar(im, ax = ax)

                hfield = np.sqrt(np.abs(hfield))

                if run < 0:
                    h_ref.append(np.sum(np.abs(hfield)*mag_dx*mag_dy))
                    ts.append(snap*0.5)
                else:
                    h_all[ri-1].append(np.sum(np.abs(hfield)*mag_dx*mag_dy))

        plt.tight_layout()
        plt.close()

        np.save('constant_omega/h_all.npy', np.array(h_all))
        np.save('constant_omega/h_ref.npy', np.array(h_ref))
        np.save('constant_omega/ts.npy', np.array(ts))
        np.save('constant_omega/omegas.npy', np.array(omegas))

nruns = 7

if os.uname()[1] == 'brillouin.dur.ac.uk':
    os.system('scp -r trcn27@hamilton8.dur.ac.uk:/home/trcn27/mf_emergence/writeup/constant_omega/')
    pass

else:
    compute_inplane_helicity(0,nruns)

#compute_inplane_helicity(0,1)

#Copy data, load and plot

h_all = np.load('constant_omega/h_all.npy')
h_ref = np.load('constant_omega/h_ref.npy')
ts = np.load('constant_omega/ts.npy')
omegas = np.load('constant_omega/omegas.npy')

fig = plt.figure(figsize = (10,7))
for ri in range(nruns):
    plt.plot(ts[:len(h_all[ri])], h_all[ri], label = ('Omega factor= %.4f' % omegas[ri]))

plt.plot(ts, h_ref, c= 'black', linestyle = 'dashed', label = 'LARE Reference')

plt.legend()

plt.xlabel('Time')
plt.ylabel('In-plane helicity A.B')
plt.tight_layout()
plt.savefig('constant_omega.png')
plt.savefig('constant_omega.pdf')

plt.show()
