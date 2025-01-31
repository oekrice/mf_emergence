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
from numpy import random
import pyvista as pv


matplotlib.rcParams['text.usetex'] = True

cs = plt.cm.plasma(np.linspace(0.1,0.9,7))
fig_width = 1.0*(513.11743/72)

#There are seven runs. Use the same script here and on Hamilton for neatness purposes.

class trace_fieldlines():
    def __init__(self, run, snap_min, snap_max, show = 0):

        data_directory = ('./constant_omega/%03d/' % run)
        #Establish random seeds for field line plotting
        if not os.path.isfile('start_seeds.txt'):
            self.start_seeds = random.rand(1000**2)
            np.savetxt('start_seeds.txt', self.start_seeds, delimiter = ',')
        else:
            self.start_seeds = np.loadtxt('start_seeds.txt', delimiter = ',')

        self.show = show
        self.machine_flag = 0
        print('Data root', data_directory)
        #Establish grid parameters (can be read in from elsewhere of course)
        for snap_number in range(snap_min, snap_max):
            self.run = run
            self.snap = snap_number
            self.print_flag = 1

            self.save_number = self.snap
            self.option = 2   #tracing a plotting options (1 for jet, 2 for emergence)

            self.data_source = 0.
            self.data_root = data_directory
            data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
            self.bx = np.swapaxes(data.variables['bx'][:],0,2)
            self.by = np.swapaxes(data.variables['by'][:],0,2)
            self.bz = np.swapaxes(data.variables['bz'][:],0,2)
            data.close()

            #Establish start points for the field line plotting
            self.max_line_length = 10000
            self.ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
            self.weakness_limit = 1e-4  #Minimum field strength to stop plotting
            self.line_plot_length = 100  #To save time while plotting, reduce the length of the plotted lines

            #Import bz as a test of the resolutions (and for the pyvista plot)
            self.nx = np.shape(self.bz)[0]
            self.ny = np.shape(self.bz)[1]
            self.nz = np.shape(self.bz)[2] - 1

            self.x0 = -130.; self.x1 = 130.
            self.y0 = -130.; self.y1 = 130.
            self.z0 = 0.; self.z1 = 130.

            self.xs = np.linspace(self.x0,self.x1,self.nx+1)
            self.ys = np.linspace(self.y0,self.y1,self.ny+1)
            self.zs = np.linspace(self.z0,self.z1,self.nz+1)

            self.xc = np.zeros(self.nx + 2)
            self.yc = np.zeros(self.ny + 2)
            self.zc = np.zeros(self.nz + 2)

            self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
            self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
            self.zc[1:-1] = 0.5*(self.zs[1:] + self.zs[:-1])

            self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
            self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])
            self.zc[0] = self.zc[0] - (self.zc[2] - self.zc[2])

            self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
            self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])
            self.zc[-1] = self.zc[-2] + (self.zc[-2] - self.zc[-3])

            self.dx = self.xs[1] - self.xs[0]
            self.dy = self.ys[1] - self.ys[0]
            self.dz = self.zs[1] - self.zs[0]

            #Folder admin
            if not os.path.exists('./fl_data/'):
                os.mkdir('fl_data')

            #flh = FLH(self)    #Do field-line helicity things
            #self.flh_photo = flh.flh_photo #FLH density on the photosphere
            #self.plot_base = self.flh_photo
            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            self.plot_base = self.bz[:,:,0]
            #self.plot_base = np.ones((self.nx, self.ny))
            #Find start points
            self.set_starts()

            #Create runtime variables for fortran
            self.setup_tracer()
            #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
            self.trace_lines_fortran()
            #Plot the field lines (using pyvista)
            if not os.path.exists('./plots/'):
                os.mkdir('plots')

            os.system('rm ./fl_data/flines%03d.nc' % self.snap)
            os.system('rm ./fl_data/flparameters%03d.txt' % self.snap)
            os.system('rm ./fl_data/starts%03d.txt' % self.snap)

            self.plot_vista()


    def plot_vista(self):
        print('Plotting...')

        x, y = np.meshgrid(self.xs, self.ys)
        z = 10*np.ones((np.shape(x)))
        surface = pv.StructuredGrid(x, y, z)

        print(len(self.lines), ' lines imported')
        for li, line in enumerate(self.lines):
            line_length = 0
            for k in range(len(line)):
                if line[k,2] < 1e6 and line[k,2] >= 10.0:
                    line_length += 1
                else:
                    break
            line = line[:line_length,:]
            #Thin out the lines (if required)
            if line_length > 1:
                thin_fact = max(int(line_length/self.line_plot_length), 1)
                thinned_line = line[:line_length:thin_fact].copy()
                thinned_line[-1] = line[line_length-1].copy()
            else:
                continue

            line = np.array(thinned_line)
            doplot = True
            if line_length == 0:
                doplot = False
            #Filter if out of a specific range
            maxs = np.max(np.sqrt(line[:,0]**2 + line[:,1]**2))

            if maxs > 0.6*self.xs[-1]:
                doplot = False
            #if line[-1][2] > 11.0 and line[-1][2] < 99.0:
            #    doplot = False
            line = line.tolist()
            if doplot:
                p.add_mesh(pv.Spline(line, len(line)),color='white',line_width=0.25)

        p.add_mesh(surface, scalars= np.abs(self.plot_base), show_edges=False,cmap = 'plasma')
        p.camera.position = (200.0,100,150.0)
        p.camera.focal_point = (0,0,50)
        p.remove_scalar_bar()

    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.

        #Plot up from the surface, based on some threshold of how strong the magnetic field is... Could be fun.
        z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

        #surface_array = self.bz[:,:,z_photo]   #Distribution of surface magnetic field. This array is all that needs changing
        surface_array = self.plot_base   #Distribution of surface FLH
        #Want to plot the lines based on where the flh differences are highest

        max_surface = np.max(np.abs(surface_array)) + 1e-6

        nlines = 250

        alpha = 2.0
        alphasum = np.sum(np.abs(surface_array)**alpha)
        pb = max_surface**alpha*nlines/alphasum

        print('prob', pb, nlines, alphasum)

        self.starts = []

        cellcount = 0
        for i in range(self.nx):  #run through grid cells
            for j in range(self.ny):
                prop = np.abs(surface_array[i,j])/max_surface
                if self.start_seeds[cellcount] < pb*prop**alpha:
                    self.starts.append([self.xc[i+1],self.yc[j+1],10.0])
                cellcount += 1

        print('Tracing', len(self.starts), 'lines')

        self.nstarts = len(self.starts)
        self.starts = np.array(self.starts).reshape(self.nstarts*3)

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.05 #Tracing 'timestep' as a proportion of the grid size
        print_flag = 1  #Print some things as the tracing happens

        variables = np.zeros((30))

        variables[0] = self.run
        variables[1] = self.nx
        variables[2] = self.ny
        variables[3] = self.nz
        variables[4] = self.x0
        variables[5] = self.x1
        variables[6] = self.y0
        variables[7] = self.y1
        variables[8] = self.z0
        variables[9] = self.z1
        variables[10] = self.snap
        variables[11] = self.nstarts
        variables[12] = self.print_flag
        variables[13] = self.max_line_length
        variables[14] = self.ds
        variables[15] = self.weakness_limit
        variables[16] = self.data_source
        variables[17] = self.machine_flag

        np.savetxt('./fl_data/flparameters%03d.txt' % self.snap, variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('./fl_data/starts%03d.txt' % self.snap, self.starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):
        os.system('make')
        if os.uname()[1] == 'brillouin.dur.ac.uk' or os.uname()[1] == 'modigliani.dur.ac.uk':
            os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
            os.system('mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        else:
            os.system('mpirun -np 1 ./bin/fltrace %d' % self.snap)

        try:
            data = netcdf_file('./fl_data/flines%03d.nc' % self.snap, 'r', mmap=False)
            print('Field lines found')

        except:
            print('File not found')

        self.lines = np.swapaxes(data.variables['lines'][:],0,2)

def niceformat(number):
    #Gives a string that looks like 5 x 10^-4 or whatever, based on the significant figures of the input

    mag = 10   #Magnitude to care about
    while mag > -16:
        if int(number / 10**mag) > 0:
            #This is the correct magnitude
            main = int(number / 10**mag)
            exp = mag
            break
        else:
            mag -= 1
    output = str(main) + 'e' + str(mag)
    return output


nruns = 7
snap_num = 33

runs = [0,2,4,6]#,3]#,5,7]

omegas = np.load('constant_omega/omegas.npy')

if os.uname()[1] == 'brillouin.dur.ac.uk':
    for run in runs:
        try:
            os.mkdir('constant_omega/%03d/' % run)
        except:
            pass
        #Copy files across from Hamilton
        #os.system('scp -r trcn27@hamilton8.dur.ac.uk:/nobackup/trcn27/mf3d0/%03d/%04d.nc ./constant_omega/%03d/' % (run, snap_num, run))


    show = False
    p = pv.Plotter(off_screen=not(show), shape = (2,2))
    p.background_color = "black"

    for ri, run in enumerate(runs):
        p.subplot(ri//2, ri%2)

        trace_fieldlines(run, snap_num, snap_num + 1, show = 1)  #This contains the main plotting bit but NOT the preamble or postamble

        #p.add_title('Omega = %s' % niceformat(omegas[run]), color = 'Red', font_size = 12)

        p.add_text('Omega = %s' % niceformat(omegas[run]), position = 'lower_edge', color = 'Green', font_size = 12)

    if show:
        p.show()
    else:
        p.show(screenshot='fltrace_constant.png', window_size = (1000,1000))
