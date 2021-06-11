#!/usr/bin/python
# -*- coding: utf-8 -*-

import h5py
import matplotlib
import numpy as np
import pylab as P
from mpl_toolkits.axes_grid1 import AxesGrid
import plot_utils as pu
import read_dataset as rd
matplotlib.use('cairo')      # choose output format

def draw_particles(ax, p1, p2, pm, nbins, min1, max1, min2, max2, drawd, pcolor):
    if nbins > 1:
        ax.hist2d(p1, p2, nbins, weights=pm, range=[[min1, max1], [min2, max2]], norm=matplotlib.colors.LogNorm(), cmap=pcolor)
        if not drawd:
            ax.set_facecolor('xkcd:black')
    else:
        ax.scatter(p1, p2, c=pcolor, marker=".")
        ax.set_xlim(min1, max1)
        ax.set_ylim(min2, max2)
    return ax


def plotcompose(pthfilen, var, output, options):
    umin, umax, cmap, pcolor, sctype, cu, cx, cy, cz, drawd, drawp, nbins, uaxes = options
    h5f = h5py.File(pthfilen, 'r')
    time = h5f.attrs['time'][0]
    utim = h5f['dataset_units']['time_unit'].attrs['unit']
    ulenf = h5f['dataset_units']['length_unit'].attrs['unit']
    usc, ulen, uupd = rd.change_units(ulenf, uaxes)
    if drawd:
        uvar = h5f['dataset_units'][var].attrs['unit']
    nx, ny, nz = h5f['simulation_parameters'].attrs['domain_dimensions']
    xmin, ymin, zmin = h5f['simulation_parameters'].attrs['domain_left_edge']
    xmax, ymax, zmax = h5f['simulation_parameters'].attrs['domain_right_edge']
    if uupd:
        xmin, ymin, zmin = xmin / usc, ymin / usc, zmin / usc
        xmax, ymax, zmax = xmax / usc, ymax / usc, zmax / usc
    h5f.close()

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    if drawp:
        px, py, pz, pm = rd.collect_particles(pthfilen, nbins)
        if uupd:
            px, py, pz = px / usc, py / usc, pz / usc

    if drawd:
        dset = rd.collect_dataset(pthfilen, var)

        ix = int(nx / 2 - 1)
        iy = int(ny / 2 - 1)
        iz = int(nz / 2 - 1)
        if cu:
            ix = int(min(nx - 1, max(0, np.floor(nx * (cx - xmin) / (xmax - xmin)))))
            iy = int(min(ny - 1, max(0, np.floor(ny * (cy - ymin) / (ymax - ymin)))))
            iz = int(min(nz - 1, max(0, np.floor(nz * (cz - zmin) / (zmax - zmin)))))
            print('Ordered plot center', cx, cy, cz, ' gives following uniform grid indices:', ix, iy, iz)

        xy = dset[:, :, iz]
        xz = dset[:, iy, :].swapaxes(0, 1)
        yz = dset[ix, :, :].swapaxes(0, 1)
        d3min, d3max = np.min(dset), np.max(dset)
        d2max = max(np.max(xz), np.max(xy), np.max(yz))
        d2min = min(np.min(xz), np.min(xy), np.min(yz))

        xy, xz, yz, vmin, vmax = pu.scale_manage(sctype, xy, xz, yz, umin, umax, d2min, d2max)

        print('3D data value range: ', d3min, d3max)
        print('Slices  value range: ', d2min, d2max)
        print('Plotted value range: ', vmin, vmax)
    fig = P.figure(1, figsize=(10, 10.5))

    if drawd:
        cbar_mode = 'single'
    else:
        cbar_mode = 'none'

    grid = AxesGrid(fig, 111, nrows_ncols=(2, 2), axes_pad=0.2, aspect=True, cbar_mode=cbar_mode, label_mode="L",)

    ax = grid[3]
    if drawd:
        a = ax.imshow(xz, origin="lower", extent=[xmin, xmax, zmin, zmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax = draw_particles(ax, px, pz, pm, nbins, xmin, xmax, zmin, zmax, drawd, pcolor)
    ax.set_xlabel("x [%s]" % pu.labelx()(ulen))
    ax.set_ylabel("z [%s]" % pu.labelx()(ulen))

    ax = grid[0]
    if drawd:
        a = ax.imshow(xy, origin="lower", extent=[ymin, ymax, xmin, xmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax = draw_particles(ax, py, px, pm, nbins, ymin, ymax, xmin, xmax, drawd, pcolor)
    ax.set_ylabel("x [%s]" % pu.labelx()(ulen))
    ax.set_xlabel("y [%s]" % pu.labelx()(ulen))
    ax.set_title(timep)

    ax = grid[2]
    if drawd:
        a = ax.imshow(yz, origin="lower", extent=[ymin, ymax, zmin, zmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax = draw_particles(ax, py, pz, pm, nbins, ymin, ymax, zmin, zmax, drawd, pcolor)
    ax.set_xlabel("y [%s]" % pu.labelx()(ulen))
    ax.set_ylabel("z [%s]" % pu.labelx()(ulen))

    if drawd:
        bar = grid.cbar_axes[0]
        bar.axis["right"].toggle(all=True)
        cbar = P.colorbar(a, cax=bar, format='%.1e', drawedges=False)
        cbar.ax.set_ylabel(var + " [%s]" % pu.labelx()(uvar))

    P.draw()
    P.savefig(output, facecolor='white')
    print(output, "written to disk")
    P.clf()
