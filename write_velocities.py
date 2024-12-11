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

smooth = "12"
inputDir = "magnetograms"

for snap in range(0,500):

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

    #Data loaded in at the original resolution. Even if this is the same as the intended, need to interpolate as it will be one cell off (bugger)...



    fid = netcdf_file(vfield_filename, 'w')
    fid.createDimension('xc', nx)
    fid.createDimension('yc', ny)

    vid = fid.createVariable('xc', 'd', ('xc',))
    vid[:] = xc
    vid = fid.createVariable('yc', 'd', ('yc',))
    vid[:] = yc

    #Transposes are necessary as it's easier to flip here than in Fortran
    vid = fid.createVariable('vx', 'd', ('yc','xc'))
    vid[:] = dataUx
    vid = fid.createVariable('vy', 'd', ('yc','xc'))
    vid[:] = dataUy
    vid = fid.createVariable('vz', 'd', ('yc','xc'))
    vid[:] = dataUz

    fig, axs = plt.subplots(1,3)
    axs[0].pcolormesh(dataUx.T)
    axs[1].pcolormesh(dataUy.T)
    axs[2].pcolormesh(dataUz.T)

    plt.savefig('./velocity_plots/%03d.png' % snap)

    fid.close()



