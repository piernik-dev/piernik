#!/usr/bin/python
# -*- coding: utf-8 -*-

import h5py
import matplotlib
import numpy as np
import pylab as P
from mpl_toolkits.axes_grid1 import AxesGrid
import plot_utils as pu
import read_dataset as rd
import time as timer
matplotlib.use('cairo')      # choose output format


def plot_axes(ax, ulen, l1, min1, max1, l2, min2, max2):
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_xlabel("%s [%s]" % (l1, pu.labelx()(ulen)))
    ax.set_ylabel("%s [%s]" % (l2, pu.labelx()(ulen)))
    return ax


def draw_plotcomponent(ax, refis, drawd, smin, smax, vmin, vmax, cmap, ncut, n1, n2):
    if drawd:
        for blks in refis:
            for bl in blks:
                binb, bxyz, ble, bre = bl
                if binb[ncut]:
                    ag = ax.imshow(bxyz[ncut], origin="lower", extent=[ble[n1], bre[n1], ble[n2], bre[n2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    return ax, ag

def draw_particles(ax, p1, p2, pm, nbins, ranges, drawd, pcolor, psize):
    if nbins > 1:
        ah = ax.hist2d(p1, p2, nbins, weights=pm, range=[ranges[0:2], ranges[2:4]], norm=matplotlib.colors.LogNorm(), cmap=pcolor)
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
    usc, ulen, uupd = pu.change_units(ulenf, uaxes)
    if drawd:
        uvar = h5f['dataset_units'][var].attrs['unit']
    if drawh:
        umass = h5f['dataset_units']['mass_unit'].attrs['unit']
    nd = h5f['simulation_parameters'].attrs['domain_dimensions']
    smin = h5f['simulation_parameters'].attrs['domain_left_edge']
    smax = h5f['simulation_parameters'].attrs['domain_right_edge']
    if uupd:
        smin = pu.list3_division(smin, usc)
        smax = pu.list3_division(smax, usc)
    cgcount = int(h5f['data'].attrs['cg_count'])
    glevels = h5f['grid_level'][:]
    maxglev = max(glevels)

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    if drawp:
        px, py, pz, pm = rd.collect_particles(h5f, nbins)
        if uupd:
            px, py, pz = px / usc, py / usc, pz / usc

    if drawd:
        if not cu:
            center = (smax[0] + smin[0]) / 2.0, (smax[1] + smin[1]) / 2.0, (smax[2] + smin[2]) / 2.0

        xy, xz, yz, extr = rd.reconstruct_uniform(h5f, var, cu, center, nd, smin, smax)
        d2min, d2max, d3min, d3max = extr
        block = [True, True, True], [yz, xz, xy], smin, smax
        refis = [[block, ], ]

        refis = rd.collect_gridlevels(h5f, var, maxglev, cgcount, center, usc)

        xy, xz, yz, vmin, vmax = pu.scale_manage(sctype, xy, xz, yz, umin, umax, d2min, d2max)

        print('3D data value range: ', d3min, d3max)
        print('Slices  value range: ', d2min, d2max)
        print('Plotted value range: ', vmin, vmax)

    h5f.close()

    fig = P.figure(1, figsize=(10, 10.5))

    cbar_mode = pu.colorbar_mode(drawd, drawh)

    if not zoom[0]:
        zoom = False, smin, smax

    grid = AxesGrid(fig, 111, nrows_ncols=(2, 2), axes_pad=0.2, aspect=True, cbar_mode=cbar_mode, label_mode="L",)

    ax = grid[3]
    ax, ag = draw_plotcomponent(ax, refis, drawd, smin, smax, vmin, vmax, cmap, 1, 0, 2)
    if drawp:
        extent = [smin[0], smax[0], smin[2], smax[2]]
        ax, ah = draw_particles(ax, px, pz, pm, nbins, extent, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "x", zoom[1][0], zoom[2][0], "z", zoom[1][2], zoom[2][2])

    ax = grid[0]
    ax, ag = draw_plotcomponent(ax, refis, drawd, smin, smax, vmin, vmax, cmap, 2, 1, 0)
    if drawp:
        extent = [smin[1], smax[1], smin[0], smax[0]]
        ax, ah = draw_particles(ax, py, px, pm, nbins, extent, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "y", zoom[1][1], zoom[2][1], "x", zoom[1][0], zoom[2][0])
    ax.set_title(timep)

    ax = grid[2]
    ax, ag = draw_plotcomponent(ax, refis, drawd, smin, smax, vmin, vmax, cmap, 0, 1, 2)
    if drawp:
        extent = [smin[1], smax[1], smin[2], smax[2]]
        ax, ah = draw_particles(ax, py, pz, pm, nbins, extent, drawd, pcolor, psize)
    ax = plot_axes(ax, ulen, "y", zoom[1][1], zoom[2][1], "z", zoom[1][2], zoom[2][2])

    if drawh:
        add_cbar(cbar_mode, grid, ah[3], 0.17, 'particle mass histogram' + " [%s]" % pu.labelx()(umass))

    if drawd:
        add_cbar(cbar_mode, grid, ag, 0.23, var + " [%s]" % pu.labelx()(uvar))

    P.draw()
    P.savefig(output, facecolor='white')
    print(output, "written to disk")
    P.clf()
