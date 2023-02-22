#!/usr/bin/python
# -*- coding: utf-8 -*-

import h5py
import matplotlib
import numpy as np
import pylab as P
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plot_utils as pu
import pvf_settings as ps
import read_dataset as rd
import time as timer
# matplotlib.use('cairo')      # choose output format


def plot_axes(ax, ulen, l1, min1, max1, l2, min2, max2):
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_xlabel("%s [%s]" % (l1, pu.labelx()(ulen)))
    ax.set_ylabel("%s [%s]" % (l2, pu.labelx()(ulen)))
    return ax


def plot1d(refis, field, parts, equip1d, ncut, n1, n2):
    smin, smax, zoom, ulen, sctype, umin, umax, linstyl, output, timep = equip1d
    fig1d = P.figure(ncut + 2, figsize=(10, 8))
    ax = fig1d.add_subplot(111)
    P.xlim(zoom[1][ncut], zoom[2][ncut])
    label = []
    axis = "xyz"[ncut]
    hybrid_plot = False
    if field[0] and parts[0]:
        hybrid_plot = (parts[3] > 1)

    if field[0]:
        vmin, vmax, symmin, cmap, autsc, labf = field[1:]
        label.append(labf)
        if not autsc and not hybrid_plot:
            P.ylim(vmin, vmax)

        for blks in refis:
            for bl in blks:
                binb, ble, bre, level, b1d = bl[1:]
                if binb[n1] and binb[n2]:
                    if b1d != []:
                        if hybrid_plot:
                            bplot = b1d[ncut]
                        else:
                            bplot = pu.scale_plotarray(b1d[ncut], sctype, symmin)
                        dxh = (bre[ncut] - ble[ncut]) / float(len(b1d[ncut])) / 2.0
                        vax = np.linspace(ble[ncut] + dxh, bre[ncut] - dxh, len(b1d[ncut]))
                        ax.plot(vax, bplot, linstyl[level], color=ps.plot1d_linecolor)

    if hybrid_plot:
        hl = pu.scale_translate(sctype, vmin, vmax, symmin, True)
    if not field[0]:
        autsc = (umin == umax)
        hl = pu.scale_translate(sctype, umin, umax, np.abs(umin), False)
    if hybrid_plot or not field[0]:
        if hl[0] == 3:
            if hl[3] > 0.:
                ax.set_yscale('symlog', linthresh=hl[3])
            else:
                ax.set_yscale('symlog')
        elif hl[0] == 2:
            ax.set_yscale('log')
        if not autsc:
            P.ylim(hl[1], hl[2])

    if parts[0]:
        pxyz, pm, nbins, pcolor, psize, player, pstype, labh = parts[1:]
        pn1, pmm = pxyz[ncut], pm
        ax = plot1d_particles(ax, pn1, pmm, nbins, [smin[ncut], smax[ncut]], pcolor, psize)
        label.append(labh)

    P.ylabel('  |  '.join(label))
    P.xlabel("%s [%s]" % (axis, pu.labelx()(ulen)))
    P.title(timep)
    P.tight_layout()
    P.draw()
    out1d = output[0] + axis + '_' + output[1]
    fig1d.savefig(out1d, facecolor=ps.f_facecolor)
    print(out1d, "written to disk")
    P.clf()


def draw_plotcomponent(ax, refis, field, parts, equip2d, ncut, n1, n2):
    smin, smax, zoom, ulen, sctype, drawg, gcolor, center = equip2d
    ag, ah = [], []
    if field[0] or drawg:
        if field[0]:
            vmin, vmax, symmin, cmap = field[1:-2]
        for blks in refis:
            for bl in blks:
                bxyz, binb, ble, bre, level = bl[:-1]
                if binb[ncut]:
                    if bxyz != []:
                        bplot = pu.scale_plotarray(bxyz[ncut], sctype, symmin)
                        ag = ax.imshow(bplot, origin="lower", extent=[ble[n1], bre[n1], ble[n2], bre[n2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
                    if drawg:
                        ax.plot([ble[n1], ble[n1], bre[n1], bre[n1], ble[n1]], [ble[n2], bre[n2], bre[n2], ble[n2], ble[n2]], '-', linewidth=ps.grid_linewidth, alpha=0.1 * float(level + 1), color=gcolor[level], zorder=4)
    if parts[0]:
        pxyz, pm, nbins, pcolor, psize, player, pstype = parts[1:-1]
        if player[0] and player[ncut + 1] != '0':
            pmask = np.abs(pxyz[ncut] - center[ncut]) <= float(player[ncut + 1])
            pn1, pn2 = pxyz[n1][pmask], pxyz[n2][pmask]
            if nbins > 1:
                pmm = pm[pmask]
            else:
                pmm = pm
        else:
            pn1, pn2, pmm = pxyz[n1], pxyz[n2], pm
        ax, ah = draw_particles(ax, pn1, pn2, pmm, nbins, [smin[n1], smax[n1], smin[n2], smax[n2]], field[0], pcolor, psize, pstype)
    ax = plot_axes(ax, ulen, "xyz"[n1], zoom[1][n1], zoom[2][n1], "xyz"[n2], zoom[1][n2], zoom[2][n2])
    ax.set_xticks([center[n1]], minor=True)
    ax.set_yticks([center[n2]], minor=True)
    ax.tick_params(axis='x', which='minor', color=ps.centertick_color, bottom='on', top='on', width=ps.centertick_width, length=ps.centertick_length)
    ax.tick_params(axis='y', which='minor', color=ps.centertick_color, left='on', right='on', width=ps.centertick_width, length=ps.centertick_length)
    return ax, ag, ah


def draw_particles(ax, p1, p2, pm, nbins, ranges, drawd, pcolor, psize, pstype):
    if nbins > 1:
        if pstype[0]:
            norm = matplotlib.colors.LogNorm(vmin=pstype[1], vmax=pstype[2])
        else:
            norm = matplotlib.colors.Normalize(vmin=pstype[1], vmax=pstype[2])
        ah = ax.hist2d(p1, p2, nbins, weights=pm, range=[ranges[0:2], ranges[2:4]], norm=norm, cmap=pcolor, alpha=ps.hist2d_alpha)
        if not drawd:
            ax.set_facecolor(ps.hist2d_facecolor)
    else:
        ah = []
        if psize <= 0:
            psize = matplotlib.rcParams['lines.markersize']**2
        ax.scatter(p1, p2, c=pcolor, marker=ps.particles_marker, s=psize)
    return ax, ah


def plot1d_particles(ax, p1, pm, nbins, ranges, pcolor, psize):
    if nbins > 1:
        ax.hist(p1, nbins, weights=pm, range=ranges, color=ps.particles_color, alpha=ps.hist2d_alpha)
    else:
        if psize <= 0:
            psize = matplotlib.rcParams['lines.markersize']**2
        ax.scatter(p1, p1 * 0., c=pcolor, marker=ps.particles_marker, s=psize)
    return ax


def add_cbar(figmode, cbar_mode, grid, ab, ic, clab, sct, field):
    icb = pu.cbsplace[figmode][ic]
    clf = [ps.cbar_hist2d_label_format, ps.cbar_plot2d_label_format][ic]
    if cbar_mode == 'none':
        axg = grid[icb]
        if figmode == 4 or figmode == 5:
            axg.set_xlim((grid[icb - 1].get_xlim()[0] / 2., grid[icb - 1].get_xlim()[1] / 2.))
        pu.color_axes(axg, 'white')
        fr = [ps.cbar_hist2d_shift, ps.cbar_plot2d_shift]
        bar = inset_axes(axg, width='100%', height='100%', bbox_to_anchor=(fr[ic], 0.0, 0.06, 1.0), bbox_transform=axg.transAxes, loc=2, borderpad=0)
        cbarh = P.colorbar(ab, cax=bar, format=clf, drawedges=False)
    else:
        if cbar_mode == 'each':
            bar = grid.cbar_axes[icb]
        elif cbar_mode == 'single':
            bar = grid.cbar_axes[0]
        bar.axis["right"].toggle(all=True)
        cbarh = P.colorbar(ab, cax=bar, format=clf, drawedges=False)
    cbarh.ax.set_ylabel(clab)
    if cbar_mode == 'none':
        cbarh.ax.yaxis.set_label_coords(ps.cbar_label_coords[0], ps.cbar_label_coords[1])
    if ic == 1 and pu.whether_symlog(sct):
        slticks, sltlabs = [], []
        for tick in cbarh.get_ticks():
            if tick >= field[1] and tick <= field[2]:
                slticks.append(tick)
                sltlabs.append(clf % (np.sign(tick) * 10.**(np.abs(tick)) * field[3]))
        if (matplotlib.__version__ >= '3.5.0'):
            cbarh.set_ticks(slticks, labels=sltlabs)
        else:
            cbarh.set_ticks(slticks)
            cbarh.set_ticklabels(sltlabs)


def plotcompose(pthfilen, var, output, options):
    axc, umin, umax, cmap, pcolor, player, psize, sctype, pstype, cu, center, cmpr, drawg, drawd, drawu, drawa, drawp, nbins, uaxes, zoom, plotlevels, gridlist, gcolor, linstyl = options
    labh = ps.particles_label
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
        labh = labh[:-1] + ' mass histogram' + " [%s]" % pu.labelx()(umass)
    smin = h5f['simulation_parameters'].attrs['domain_left_edge']
    smax = h5f['simulation_parameters'].attrs['domain_right_edge']
    n_d = h5f['simulation_parameters'].attrs['domain_dimensions']
    if uupd:
        smin = pu.list3_division(smin, usc)
        smax = pu.list3_division(smax, usc)
    cgcount = int(h5f['data'].attrs['cg_count'])
    glevels = h5f['grid_level'][:]
    maxglev = max(glevels)

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    parts, field = [drawp, ], [drawd, ]

    if not cu:
        center = (smax[0] + smin[0]) / 2.0, (smax[1] + smin[1]) / 2.0, (smax[2] + smin[2]) / 2.0

    drawa, drawu = pu.choose_amr_or_uniform(drawa, drawu, drawd, drawg, drawp, maxglev, gridlist)
    plotlevels = pu.check_plotlevels(plotlevels, maxglev, pthfilen, True)
    gridlist = pu.sanitize_gridlist(gridlist, cgcount)
    unavail, cmpr, drawa, drawu = rd.manage_compare(cmpr, pthfilen, h5f, var, plotlevels, gridlist, drawa, drawu)
    if len(plotlevels) == 0 or unavail:
        return

    linstyl = pu.linestyles(linstyl, maxglev, plotlevels)
    if drawg:
        gcolor = pu.reorder_gridcolorlist(gcolor, maxglev, plotlevels)

    if drawp:
        pinfile, pxyz, pm = rd.collect_particles(h5f, drawh, center, player, uupd, usc, plotlevels, gridlist)
        parts = pinfile, pxyz, pm, nbins, pcolor, psize, player, pstype, labh
        drawh = drawh and pinfile

    draw1D, draw2D, figmode = pu.check_1D2Ddefaults(axc, n_d, drawd and drawh)

    refis = []
    if drawd or drawg:
        refis, extr, center = rd.collect_gridlevels(h5f, var, cmpr, refis, maxglev, plotlevels, gridlist, cgcount, center, usc, drawd, drawu, drawa, drawg, draw1D, draw2D)

        if refis == [] or pu.list_any(extr, []):
            drawd = False
            field = [drawd, ]
        else:
            if drawd:
                vmin, vmax, symmin, autsc = pu.scale_manage(sctype, refis, umin, umax, draw1D, draw2D, extr)

                vlab = pu.labellog(sctype, symmin, cmpr[0]) + var + pu.manage_units(uvar)
                field = drawd, vmin, vmax, symmin, cmap, autsc, vlab

    zoom = rd.level_zoom(h5f, gridlist, zoom, smin, smax)

    h5f.close()
    if cmpr[0]:
        cmpr[2].close()

    if not (parts[0] or drawd or drawg):
        print('No particles or datafields or grids to plot. Skipping.')
        return

    cbar_mode = pu.colorbar_mode(drawd, drawh, figmode)

    equip1d = smin, smax, zoom, ulen, sctype, umin, umax, linstyl, output, timep
    equip2d = smin, smax, zoom, ulen, sctype, drawg, gcolor, center

    if draw1D[0]:
        plot1d(refis, field, parts, equip1d, 0, 1, 2)
    if draw1D[1]:
        plot1d(refis, field, parts, equip1d, 1, 0, 2)
    if draw1D[2]:
        plot1d(refis, field, parts, equip1d, 2, 0, 1)

    if any(draw2D):
        fig = P.figure(1, figsize=pu.figsizes[figmode])

        grid = AxesGrid(fig, 111, nrows_ncols=pu.figrwcls[figmode], axes_pad=ps.plot2d_axes_pad, aspect=True, cbar_mode=cbar_mode, label_mode="L",)
        ag0, ag2, ag3 = [], [], []

        if draw2D[0]:
            ax = grid[pu.figplace[figmode][0]]
            ax, ag3, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 0, 1, 2)

        if draw2D[1]:
            ax = grid[pu.figplace[figmode][1]]
            ax, ag2, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 1, 0, 2)

        if draw2D[2]:
            ax = grid[pu.figplace[figmode][2]]
            ax, ag0, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 2, 0, 1)
        ax.set_title(timep)

        if drawh:
            add_cbar(figmode, cbar_mode, grid, ah[3], 0, labh, sctype, field)

        if drawd:
            add_cbar(figmode, cbar_mode, grid, pu.take_nonempty([ag0, ag2, ag3]), 1, vlab, sctype, field)

        P.draw()
        out2d = output[0] + pu.plane_in_outputname(figmode, draw2D) + output[1]
        P.savefig(out2d, facecolor=ps.f_facecolor)
        print(out2d, "written to disk")
        P.clf()
