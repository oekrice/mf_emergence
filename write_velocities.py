#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Read in the 'magnetograms' from the .nc files
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.interpolate import RegularGridInterpolator

smooth = "12"
inputDir = "magnetograms"

#Define the grid onto which these velocities should be interpolated
nx_output = 64
ny_output = 64

for snap in range(0,499):

    print('Writing snap', snap)
    vfield_filename = './magnetograms/velocity%04d.nc' % snap

    dataUx = np.loadtxt(inputDir+"/Ux_"+smooth+"_"+str(snap)+".txt")

    nx_import = int(np.sqrt(np.size(dataUx)))
    ny_import = int(np.sqrt(np.size(dataUx)))

    xs = np.linspace(0,1,nx_import+1)
    ys = np.linspace(0,1,ny_import+1)

    xc = 0.5*(xs[1:] + xs[:-1])
    yc = 0.5*(ys[1:] + ys[:-1])

    dataUx = np.loadtxt(inputDir+"/Ux_"+smooth+"_"+str(snap)+".txt")
    dataUy = np.loadtxt(inputDir+"/Uy_"+smooth+"_"+str(snap)+".txt")
    dataUz = np.loadtxt(inputDir+"/Uz_"+smooth+"_"+str(snap)+".txt")

    dataUx = np.reshape(dataUx,(nx_import,ny_import))
    dataUy = np.reshape(dataUy,(nx_import,ny_import))
    dataUz = np.reshape(dataUz,(nx_import,ny_import))

    dataBx = np.loadtxt(inputDir+"/bx_"+str(snap)+".txt")
    dataBy = np.loadtxt(inputDir+"/by_"+str(snap)+".txt")
    dataBz = np.loadtxt(inputDir+"/bz_"+str(snap)+".txt")

    dataBx = np.reshape(dataBx,(nx_import,ny_import))
    dataBy = np.reshape(dataBy,(nx_import,ny_import))
    dataBz = np.reshape(dataBz,(nx_import,ny_import))


    #Data loaded in at the original resolution. Even if this is the same as the intended, need to interpolate as it will be one cell off (bugger)...

    xs_output = np.linspace(0,1,nx_output+1)
    ys_output = np.linspace(0,1,ny_output+1)

    xc_output = 0.5*(xs_output[1:] + xs_output[:-1])
    yc_output = 0.5*(ys_output[1:] + ys_output[:-1])

    X, Y = np.meshgrid(xs_output, ys_output)

    interp_x = RegularGridInterpolator((xc, yc), dataUx, bounds_error = False, method = 'linear', fill_value = None)
    interp_y = RegularGridInterpolator((xc, yc), dataUy, bounds_error = False, method = 'linear', fill_value = None)
    interp_z = RegularGridInterpolator((xc, yc), dataUz, bounds_error = False, method = 'linear', fill_value = None)

    vx_out = interp_x((X, Y))
    vy_out = interp_y((X, Y))
    vz_out = interp_z((X, Y))

    interp_x = RegularGridInterpolator((xc, yc), dataBx, bounds_error = False, method = 'linear', fill_value = None)
    interp_y = RegularGridInterpolator((xc, yc), dataBy, bounds_error = False, method = 'linear', fill_value = None)
    interp_z = RegularGridInterpolator((xc, yc), dataBz, bounds_error = False, method = 'linear', fill_value = None)

    X, Y = np.meshgrid(xs_output, yc_output)
    bx_out = interp_x((X, Y))
    X, Y = np.meshgrid(xc_output, ys_output)
    by_out = interp_y((X, Y))
    X, Y = np.meshgrid(xc_output, yc_output)
    bz_out = interp_z((X, Y))

    fid = netcdf_file(vfield_filename, 'w')
    fid.createDimension('xs', nx_output+1)
    fid.createDimension('ys', ny_output+1)
    fid.createDimension('xc', nx_output)
    fid.createDimension('yc', ny_output)

    vid = fid.createVariable('xs', 'd', ('xs',))
    vid[:] = xs_output
    vid = fid.createVariable('ys', 'd', ('ys',))
    vid[:] = ys_output

    vid = fid.createVariable('xc', 'd', ('xc',))
    vid[:] = xc_output
    vid = fid.createVariable('yc', 'd', ('yc',))
    vid[:] = yc_output

    #Transposes are necessary as it's easier to flip here than in Fortran
    vid = fid.createVariable('vx', 'd', ('ys','xs'))
    vid[:] = vx_out
    vid = fid.createVariable('vy', 'd', ('ys','xs'))
    vid[:] = vy_out
    vid = fid.createVariable('vz', 'd', ('ys','xs'))
    vid[:] = vz_out

    #Transposes are necessary as it's easier to flip here than in Fortran
    vid = fid.createVariable('bx', 'd', ('yc','xs'))
    vid[:] = bx_out
    vid = fid.createVariable('by', 'd', ('ys','xc'))
    vid[:] = by_out
    vid = fid.createVariable('bz', 'd', ('yc','xc'))
    vid[:] = bz_out

    print('Maxs', np.max(np.abs(vx_out)), np.max(np.abs(vy_out)), np.max(np.abs(vz_out)))
    fig, axs = plt.subplots(2,3)
    axs[0,0].pcolormesh(vx_out.T)
    axs[0,1].pcolormesh(vy_out.T)
    axs[0,2].pcolormesh(vz_out.T)

    axs[1,0].pcolormesh(bx_out.T)
    axs[1,1].pcolormesh(by_out.T)
    axs[1,2].pcolormesh(bz_out.T)


    plt.savefig('./velocity_plots/%03d.png' % snap)

    plt.close()
    fid.close()



