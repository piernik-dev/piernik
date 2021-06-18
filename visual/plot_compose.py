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


def plot_axes(ax, ulen, l1, min1, max1, l2, min2, max2):
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_xlabel("%s [%s]" % (l1, pu.labelx()(ulen)))
    ax.set_ylabel("%s [%s]" % (l2, pu.labelx()(ulen)))
    return ax


def draw_particles(ax, p1, p2, pm, nbins, min1, max1, min2, max2, drawd, pcolor, psize):
    if nbins > 1:
        ah = ax.hist2d(p1, p2, nbins, weights=pm, range=[[min1, max1], [min2, max2]], norm=matplotlib.colors.LogNorm(), cmap=pcolor)
        if not drawd:
            ax.set_facecolor('xkcd:black')
    else:
        ah = []
        if psize <= 0:
            psize = matplotlib.rcParams['lines.markersize']**2
        ax.scatter(p1, p2, c=pcolor, marker=".", s=psize)
    return ax, ah


def add_cbar(cbar_mode, grid, ab, fr, clab):
    if cbar_mode == 'none':
        bar = grid[1]
        pu.color_axes(bar, 'white')
        cbarh = P.colorbar(ab, ax=bar, format='%.1e', drawedges=False, shrink=0.4668, fraction=fr, anchor=(0.0, 0.961))
    else:
        bar = grid.cbar_axes[0]
        bar.axis["right"].toggle(all=True)
        cbarh = P.colorbar(ab, cax=bar, format='%.1e', drawedges=False)
        cbarh = P.colorbar(ab, cax=bar, format='%.1e', drawedges=False)
    cbarh.ax.set_ylabel(clab)
    if cbar_mode == 'none':
        cbarh.ax.yaxis.set_label_coords(-1.5, 0.5)


def plotcompose(pthfilen, var, output, options):
    umin, umax, cmap, pcolor, psize, sctype, cu, center, drawd, drawp, nbins, uaxes, zoom = options
    drawh = drawp and nbins > 1
    h5f = h5py.File(pthfilen, 'r')
    time = h5f.attrs['time'][0]
    utim = h5f['dataset_units']['time_unit'].attrs['unit']
    ulenf = h5f['dataset_units']['length_unit'].attrs['unit']
    usc, ulen, uupd = rd.change_units(ulenf, uaxes)
    if drawd:
        uvar = h5f['dataset_units'][var].attrs['unit']
    if drawh:
        umass = h5f['dataset_units']['mass_unit'].attrs['unit']
    nd = h5f['simulation_parameters'].attrs['domain_dimensions']
    xmin, ymin, zmin = h5f['simulation_parameters'].attrs['domain_left_edge']
    xmax, ymax, zmax = h5f['simulation_parameters'].attrs['domain_right_edge']
    if uupd:
        xmin, ymin, zmin = xmin / usc, ymin / usc, zmin / usc
        xmax, ymax, zmax = xmax / usc, ymax / usc, zmax / usc
    cgcount = int(h5f['data'].attrs['cg_count'])
    glevels = h5f['grid_level'][:]
    maxglev = max(glevels)
    h5f.close()

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    if drawp:
        px, py, pz, pm = rd.collect_particles(pthfilen, nbins)
        if uupd:
            px, py, pz = px / usc, py / usc, pz / usc

    if drawd:
        if not cu:
            center = (xmax + xmin) / 2.0, (ymax + ymin) / 2.0, (zmax + zmin) / 2.0

        xy, xz, yz, extr = rd.reconstruct_uniform(pthfilen, var, cu, center, nd, [xmin, ymin, zmin], [xmax, ymax, zmax])
        d2min, d2max, d3min, d3max = extr

        refis = rd.collect_gridlevels(pthfilen, var, maxglev, cgcount, center, usc)

        xy, xz, yz, vmin, vmax = pu.scale_manage(sctype, xy, xz, yz, umin, umax, d2min, d2max)

        print('3D data value range: ', d3min, d3max)
        print('Slices  value range: ', d2min, d2max)
        print('Plotted value range: ', vmin, vmax)

    fig = P.figure(1, figsize=(10, 10.5))

    if drawd and drawh:
        cbar_mode = 'none'
    elif drawd or drawh:
        cbar_mode = 'single'
    else:
        cbar_mode = 'none'

    if not zoom[0]:
        zoom = False, xmin, xmax, ymin, ymax, zmin, zmax

    grid = AxesGrid(fig, 111, nrows_ncols=(2, 2), axes_pad=0.2, aspect=True, cbar_mode=cbar_mode, label_mode="L",)

    ax = grid[3]
    if drawd:
        a = ax.imshow(xz, origin="lower", extent=[xmin, xmax, zmin, zmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
        for blks in refis:
            for bl in blks:
                binb, bxy, bxz, byz, ble, bre = bl
                if binb[1]:
                    a = ax.imshow(bxz, origin="lower", extent=[ble[0], bre[0], ble[2], bre[2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax, ah = draw_particles(ax, px, pz, pm, nbins, xmin, xmax, zmin, zmax, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "x", zoom[1], zoom[2], "z", zoom[5], zoom[6])

    ax = grid[0]
    if drawd:
        a = ax.imshow(xy, origin="lower", extent=[ymin, ymax, xmin, xmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
        for blks in refis:
            for bl in blks:
                binb, bxy, bxz, byz, ble, bre = bl
                if binb[2]:
                    a = ax.imshow(bxy, origin="lower", extent=[ble[1], bre[1], ble[0], bre[0]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax, ah = draw_particles(ax, py, px, pm, nbins, ymin, ymax, xmin, xmax, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "y", zoom[3], zoom[4], "x", zoom[1], zoom[2])
    ax.set_title(timep)

    ax = grid[2]
    if drawd:
        a = ax.imshow(yz, origin="lower", extent=[ymin, ymax, zmin, zmax], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
        for blks in refis:
            for bl in blks:
                binb, bxy, bxz, byz, ble, bre = bl
                if binb[0]:
                    a = ax.imshow(byz, origin="lower", extent=[ble[1], bre[1], ble[2], bre[2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if drawp:
        ax, ah = draw_particles(ax, py, pz, pm, nbins, ymin, ymax, zmin, zmax, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "y", zoom[3], zoom[4], "z", zoom[5], zoom[6])

    if drawh:
        add_cbar(cbar_mode, grid, ah[3], 0.17, 'particle mass histogram' + " [%s]" % pu.labelx()(umass))

    if drawd:
        add_cbar(cbar_mode, grid, a, 0.23, var + " [%s]" % pu.labelx()(uvar))

    P.draw()
    P.savefig(output, facecolor='white')
    print(output, "written to disk")
    P.clf()
