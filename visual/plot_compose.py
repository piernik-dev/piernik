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


def draw_plotcomponent(ax, dgrid, parts, smin, smax, zoom, ulen, ncut, n1, n2):
    ag, ah = [], []
    if dgrid[0]:
        vmin, vmax, sctype, symmin, cmap, refis = dgrid[1:]
        for blks in refis:
            for bl in blks:
                binb, bxyz, ble, bre = bl
                if binb[ncut]:
                    bplot = pu.scale_plotarray(bxyz[ncut], sctype, symmin)
                    ag = ax.imshow(bplot, origin="lower", extent=[ble[n1], bre[n1], ble[n2], bre[n2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
    if parts[0]:
        pxyz, pm, nbins, pcolor, psize = parts[1:]
        ax, ah = draw_particles(ax, pxyz[n1], pxyz[n2], pm, nbins, [smin[n1], smax[n1], smin[n2], smax[n2]], dgrid[0], pcolor, psize)
    ax = plot_axes(ax, ulen, "xyz"[n1], zoom[1][n1], zoom[2][n1], "xyz"[n2], zoom[1][n2], zoom[2][n2])
    return ax, ag, ah


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
    umin, umax, cmap, pcolor, psize, sctype, cu, center, drawd, drawu, drawa, drawp, nbins, uaxes, zoom, plotlevels, gridlist = options
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

    parts, dgrid = [drawp, ], [drawd, ]

    if drawp:
        pxyz, pm = rd.collect_particles(h5f, nbins)
        if uupd:
            pxyz = pu.list3_division(pxyz, usc)
        parts = drawp, pxyz, pm, nbins, pcolor, psize

    if drawd:
        if not cu:
            center = (smax[0] + smin[0]) / 2.0, (smax[1] + smin[1]) / 2.0, (smax[2] + smin[2]) / 2.0

        if not (drawa or drawu):
            if maxglev == 0 and gridlist == '':
                drawu = True
            else:
                drawa = True

        if plotlevels == '':
            if drawa:
                plotlevels = range(maxglev + 1)
            else:
                plotlevels = 0,
        print('LEVELS TO plot: ', plotlevels)

        if gridlist == '':
            gridlist = range(cgcount)
        else:
            grdl = []
            for ig in gridlist:
                if ig >= 0 and ig < cgcount:
                    grdl.append(ig)
                else:
                    print('Grid block %s does not exist.' % str(ig))
            gridlist = grdl
            print('GRIDLIST: ', gridlist)

        refis = []
        extr = [], [], [], []
        if drawu:
            if len(plotlevels) > 1:
                print('For uniform grid plotting only the firs given level!')
            print('Plotting base level %s' % plotlevels[0])
            refis, extr = rd.reconstruct_uniform(h5f, var, plotlevels[0], gridlist, cu, center, smin, smax)

        if drawa:
            refis, extr = rd.collect_gridlevels(h5f, var, refis, extr, maxglev, plotlevels, gridlist, cgcount, center, usc)

        if refis == []:
            drawd = False
        else:
            d2min, d2max, d3min, d3max = min(extr[0]), max(extr[1]), min(extr[2]), max(extr[3])
            vmin, vmax, symmin = pu.scale_manage(sctype, refis, umin, umax, d2min, d2max)

            print('3D data value range: ', d3min, d3max)
            print('Slices  value range: ', d2min, d2max)
            print('Plotted value range: ', vmin, vmax)
            dgrid = drawd, vmin, vmax, sctype, symmin, cmap, refis

    h5f.close()

    if not (drawp or drawd):
        print('No particles or levels to plot. Exiting.')
        return

    fig = P.figure(1, figsize=(10, 10.5))

    cbar_mode = pu.colorbar_mode(drawd, drawh)

    if not zoom[0]:
        zoom = False, smin, smax

    grid = AxesGrid(fig, 111, nrows_ncols=(2, 2), axes_pad=0.2, aspect=True, cbar_mode=cbar_mode, label_mode="L",)

    ax = grid[3]
    ax, ag, ah = draw_plotcomponent(ax, dgrid, parts, smin, smax, zoom, ulen, 1, 0, 2)

    ax = grid[0]
    ax, ag, ah = draw_plotcomponent(ax, dgrid, parts, smin, smax, zoom, ulen, 2, 1, 0)
    ax.set_title(timep)

    ax = grid[2]
    ax, ag, ah = draw_plotcomponent(ax, dgrid, parts, smin, smax, zoom, ulen, 0, 1, 2)

    if drawh:
        add_cbar(cbar_mode, grid, ah[3], 0.17, 'particle mass histogram' + " [%s]" % pu.labelx()(umass))

    if drawd:
        add_cbar(cbar_mode, grid, ag, 0.23, var + " [%s]" % pu.labelx()(uvar))

    P.draw()
    P.savefig(output, facecolor='white')
    print(output, "written to disk")
    P.clf()
