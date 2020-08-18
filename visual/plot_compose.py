#!/usr/bin/python
# -*- coding: utf-8 -*-

import h5py, matplotlib, sys
matplotlib.use('cairo')      # choose output format
import numpy as np
import os
import pylab as P
from mpl_toolkits.axes_grid1 import AxesGrid
import plot_utils as pu
import read_dataset as rd

if (len(sys.argv) < 3):
    print('PIERNIK VISUALIZATION FACILITY')
    print('Usage: ./plot_compose.py <file> <varname [varname ...]>')
    if len(sys.argv) < 2:
        exit()

pthfilen = sys.argv[1]
filen  = pthfilen.split('/')[-1]

h5f = h5py.File(pthfilen,'r')
if len(sys.argv) < 3:
    print("Available datafields: ", list(h5f['field_types'].keys()))
    exit(1)
if sys.argv[2] == "_all_":
    varlist = h5f['field_types'].keys()
else:
    varlist  = sys.argv[2:]

print(varlist)

print("Reading file: %s" % pthfilen)
for var in varlist:
    h5f = h5py.File(pthfilen,'r')
    time = h5f.attrs['time'][0]
    utim = h5f['dataset_units']['time_unit'].attrs['unit']
    ulen = h5f['dataset_units']['length_unit'].attrs['unit']
    uvar = h5f['dataset_units'][var].attrs['unit']
    nx, ny, nz = h5f['simulation_parameters'].attrs['domain_dimensions']
    xmin, ymin, zmin = h5f['simulation_parameters'].attrs['domain_left_edge']
    xmax, ymax, zmax = h5f['simulation_parameters'].attrs['domain_right_edge']
    h5f.close()

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    dset = rd.collect_dataset(pthfilen, var)

    ix = int(nx/2-1)
    iy = int(ny/2-1)
    iz = int(nz/2-1)

    xy = dset[:,:,iz]
    xz = dset[:,iy,:].swapaxes(0,1)
    yz = dset[ix,:,:].swapaxes(0,1)

    vmax = max(np.max(xz), np.max(xy), np.max(yz))
    vmin = min(np.min(xz), np.min(xy), np.min(yz))
    print('Value range: ', vmin, vmax)
    vmin, vmax = pu.fsym(vmin,vmax)

    fig = P.figure(1,figsize=(10,10.5))

    grid = AxesGrid(fig, 111, nrows_ncols = (2, 2), axes_pad = 0.2, add_all=True, aspect=True, cbar_mode='single', label_mode = "L",)

    ax = grid[3]
    a = ax.imshow(xz, origin="lower",extent=[xmin,xmax,zmin,zmax], vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_xlabel("x [%s]" % pu.labelx()(ulen))
    ax.set_ylabel("z [%s]" % pu.labelx()(ulen))

    ax = grid[0]
    a = ax.imshow(xy, origin="lower",extent=[ymin,ymax,xmin,xmax], vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_ylabel("x [%s]" % pu.labelx()(ulen))
    ax.set_xlabel("y [%s]" % pu.labelx()(ulen))
    ax.set_title(timep)

    ax = grid[2]
    a = ax.imshow(yz, origin="lower",extent=[ymin,ymax,zmin,zmax], vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_xlabel("y [%s]" % pu.labelx()(ulen))
    ax.set_ylabel("z [%s]" % pu.labelx()(ulen))

    bar = grid.cbar_axes[0]
    bar.axis["right"].toggle(all=True)
    cbar = P.colorbar(a, cax=bar,format='%.3f', drawedges=False)
    cbar.ax.set_ylabel(var+" [%s]" % pu.labelx()(uvar))

    P.draw()
    plotdir = 'frames'
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
    P.savefig(output, facecolor='white')
    print(output, "written to disk")
    P.clf()

