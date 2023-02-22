#!/usr/bin/env python
import h5py as h5
import numpy as np
import plot_utils as pu
import pvf_settings as ps


def manage_compare(cmpr, pthfilen, h5f, var, plotlevels, gridlist, drawa, drawu):
    cmpr0, cmprb, cmprf, cmprd, cmprl, cmprt, diff_struct = cmpr
    if cmpr0:
        if cmprd == '':
            cmprd = var
        if cmprl == '':
            cmprl = plotlevels
        if cmprf == '':
            cmprf = pthfilen
            if cmprl != plotlevels:
                print('Different levels to compare in the same file. Might be weird.')
            h5c = h5f
        else:
            h5c = h5.File(cmprf, 'r')

        if (cmprd not in list(h5c['field_types'].keys())):
            print(cmprd, ' is not available in the compare file ', cmprf, '. No data to compare.')
            return True, cmpr, drawa, drawu

        cmprl = pu.check_plotlevels(cmprl, max(h5c['grid_level'][:]), cmprf, False)
        diff_struct = compare_grids(h5f, h5c, plotlevels, cmprl, gridlist)
        if cmprl == []:
            return True, cmpr, drawa, drawu
        if diff_struct:
            print('Difference in datafields structure or not matching levels.')
            drawa, drawu = False, True
        cmpr = cmpr0, cmprb, h5c, cmprd, cmprl, cmprt, diff_struct
    return False, cmpr, drawa, drawu


def compare_grids(h1, h2, plotlevels, cmprl, gridlist):
    if len(plotlevels) != len(cmprl):
        return True
    cgcnt2 = cgcount = int(h2['data'].attrs['cg_count'])
    for il in range(len(plotlevels)):
        for ig in gridlist:
            h5g1 = h1['data']['grid_' + str(ig).zfill(10)]
            if ig >= cgcnt2:
                return True
            h5g2 = h2['data']['grid_' + str(ig).zfill(10)]
            if h5g1.attrs['level'] == plotlevels[il]:
                if h5g2.attrs['level'] != cmprl[il]:
                    return True
                if any(h5g1.attrs['off'] != h5g2.attrs['off']):
                    return True
                if any(h5g1.attrs['n_b'] != h5g2.attrs['n_b']):
                    return True
    return False


def reconstruct_uniform(h5f, var, cmpr, levnum, level, gridlist, center, usc, draw1D, draw2D):
    # attrs = h5f['domains']['base'].attrs
    # nd = [i * 2**level for i in attrs['n_d']]
    nd, loff, roff, ledg, redg, levelmet = frame_level(h5f, level, gridlist)
    if not levelmet:
        return False, [], []
    cmpr0, cmprb, h5c, cmprd, cmprl, cmprt, diff_struct = cmpr
    if diff_struct:
        if len(cmprl) <= levnum:
            print('Level to compare not found. Consider more detailed requirement for level comparison.')
            return False, [], []
        ndc, loc, roc, lec, rec, lmc = frame_level(h5c, cmprl[levnum], range(int(h5c['data'].attrs['cg_count'])))
        if not lmc:
            return False, [], []
        if nd == ndc and not cmprb:
            if not pu.list3_alleq(ledg, lec) or not pu.list3_alleq(redg, rec):
                print('WARNING: Edges for level %s: %s %s are different than for level %s: %s %s. Consider excluding levels or -C / --compare-adjusted-grids option.' % (level, ledg, redg, cmprl[levnum], lec, rec))
        elif cmprb:
            if (nd, ledg, redg) != (ndc, lec, rec):
                sfmin = h5f['simulation_parameters'].attrs['domain_left_edge']
                sfmax = h5f['simulation_parameters'].attrs['domain_right_edge']
                scmin = h5c['simulation_parameters'].attrs['domain_left_edge']
                scmax = h5c['simulation_parameters'].attrs['domain_right_edge']
                if not pu.list3_alleq(sfmin, scmin) or not pu.list3_alleq(sfmax, scmax):
                    print('Framing levels from domains of different edges! %s %s vs. %s %s' % (sfmin, sfmax, scmin, scmax))
                lofg = pu.list3_min(loff, loc)
                rofg = pu.list3_max(roff, roc)

                dlf = pu.list3_div(pu.list3_subtraction(redg, ledg), nd)
                dlc = pu.list3_div(pu.list3_subtraction(rec, lec), ndc)
                redg = pu.list3_add(redg, pu.list3_mult(pu.list3_subtraction(rofg, roff), dlf))
                ledg = pu.list3_add(ledg, pu.list3_mult(pu.list3_subtraction(lofg, loff), dlf))
                rec = pu.list3_add(rec, pu.list3_mult(pu.list3_subtraction(rofg, roc), dlc))
                lec = pu.list3_add(lec, pu.list3_mult(pu.list3_subtraction(lofg, loc), dlc))
                if not pu.list3_alleq(dlf, dlc):
                    print('Comparison for %s and %s may give weird results as cell sizes are different: %s %s' % (level, cmprl[levnum], dlf, dlc))

                ndgg = pu.list3_subtraction(rofg, lofg)
                nd, ndc = ndgg, ndgg
                loff, loc = lofg, lofg
        else:
            print('Comparison for levels: %s and %s not available due to unmet resolution constraints. Dimensions %s and %s do not match.' % (level, cmprl[levnum], nd, ndc))
            return False, [], []

    dset = collect_dataset(h5f, var, cmpr, level, gridlist, nd, loff)
    if diff_struct:
        dc = collect_dataset(h5c, cmprd, cmpr, cmprl[levnum], range(int(h5c['data'].attrs['cg_count'])), ndc, loc)
        if nd == ndc:
            dset = pu.execute_comparison(dset, dc, cmprt)

    if center is None:
        return pu.locate_extrema(dset, ledg, redg, nd)

    inb, ind = pu.find_indices(nd, center, ledg, redg, draw1D, draw2D, True)
    print('Plot center', center[0], center[1], center[2], 'gives indices:', ind[0], ind[1], ind[2], 'for uniform grid level', level)

    b2d, b1d, extr = take_cuts_and_lines(dset, ind, draw1D, draw2D)
    block = b2d, inb, pu.list3_division(ledg, usc), pu.list3_division(redg, usc), level, b1d

    return levelmet, block, extr


def collect_dataset(h5f, dset_name, cmpr, level, gridlist, nd, loff):
    print('Reading', dset_name)
    dset = np.full((nd[0], nd[1], nd[2]), np.nan)

    print('Reconstructing domain from cg parts')
    for ig in gridlist:
        h5g = h5f['data']['grid_' + str(ig).zfill(10)]
        if h5g.attrs['level'] == level:
            off = h5g.attrs['off'] - loff
            ngb = h5g.attrs['n_b']
            n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
            ce = n_b + off
            dset[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] = h5g[dset_name][:, :, :].swapaxes(0, 2)
            cmpr0, cmprb, h5c, cmprd, cmprl, cmprt, diff_struct = cmpr
            if cmpr0 and not diff_struct:
                dset[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] = pu.execute_comparison(dset[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]], h5c['data']['grid_' + str(ig).zfill(10)][cmprd][:, :, :].swapaxes(0, 2), cmprt)

    return dset


def frame_level(h5f, level, gridlist):
    levelmet = False
    off_started = False
    for ig in gridlist:
        h5g = h5f['data']['grid_' + str(ig).zfill(10)]
        if h5g.attrs['level'] == level:
            levelmet = True
            off = h5g.attrs['off']
            ngb = h5g.attrs['n_b']
            lft = h5g.attrs['left_edge']
            rht = h5g.attrs['right_edge']

            n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
            ce = n_b + off
            if off_started:
                lind = pu.list3_min(lind, off)
                rind = pu.list3_max(rind, ce)
                ledg = pu.list3_min(ledg, lft)
                redg = pu.list3_max(redg, rht)
            else:
                lind = off
                rind = ce
                ledg = lft
                redg = rht
                off_started = True
    nd = pu.list3_subtraction(rind, lind)
    if levelmet:
        return nd, lind, rind, ledg, redg, levelmet
    else:
        return [], [], [], [], [], False


def level_zoom(h5f, gridlist, zoom, smin, smax):
    if zoom[0]:
        if len(zoom) == 2:
            nd, loff, roff, ledg, redg, levelmet = frame_level(h5f, zoom[1], gridlist)
            if levelmet:
                zoom = True, ledg, redg
            else:
                zoom = False,
    if not zoom[0]:
        zoom = False, smin, smax
    return zoom


def collect_gridlevels(h5f, var, cmpr, refis, maxglev, plotlevels, gridlist, cgcount, center, usc, getmap, drawu, drawa, drawg, draw1D, draw2D):
    if len(center) != 3:
        curmin, curmax = np.inf, -np.inf
        locmin, locmax = None, None
        lev_num = -1
        for iref in range(maxglev + 1):
            if iref in plotlevels:
                lev_num += 1
                blks = []

                if drawu:
                    levok, bextr, cextr = reconstruct_uniform(h5f, var, cmpr, lev_num, iref, gridlist, None, usc, draw1D, draw2D)
                    if levok:
                        curmin, curmax, locmin, locmax = pu.check_extrema(curmin, curmax, locmin, locmax, bextr, cextr)

                if drawa:
                    for ib in gridlist:
                        levok, bextr, cextr = read_block(h5f, var, cmpr, ib, iref, None, usc, (getmap and drawa), draw1D, draw2D)
                        if levok:
                            curmin, curmax, locmin, locmax = pu.check_extrema(curmin, curmax, locmin, locmax, bextr, cextr)

        print('Found in position: %s %s %s \tmin value %s' % (locmin[0], locmin[1], locmin[2], curmin))
        print('Found in position: %s %s %s \tmax value %s' % (locmax[0], locmax[1], locmax[2], curmax))
        if center[0]:
            center = locmin
        elif center[1]:
            center = locmax
        print('Plot center set to: %s %s %s' % center)

    l1, h1, l2, h2, l3, h3 = [], [], [], [], [], []
    lev_num = -1
    for iref in range(maxglev + 1):
        if iref in plotlevels:
            lev_num += 1
            print('REFINEMENT ', iref)
            blks = []

            if drawu:
                levok, block, extr = reconstruct_uniform(h5f, var, cmpr, lev_num, iref, gridlist, center, usc, draw1D, draw2D)
                if levok:
                    blks.append(block)
                    if getmap:
                        l1.append(extr[0])
                        h1.append(extr[1])
                        l2.append(extr[2])
                        h2.append(extr[3])
                        l3.append(extr[4])
                        h3.append(extr[5])

            if drawa or drawg:
                for ib in gridlist:
                    levok, block, extr = read_block(h5f, var, cmpr, ib, iref, center, usc, (getmap and drawa), draw1D, draw2D)
                    if levok:
                        blks.append(block)
                        if getmap and drawa:
                            l1.append(extr[0])
                            h1.append(extr[1])
                            l2.append(extr[2])
                            h2.append(extr[3])
                            l3.append(extr[4])
                            h3.append(extr[5])
            if blks != []:
                refis.append(blks)
    return refis, [l1, h1, l2, h2, l3, h3], center


def read_block(h5f, dset_name, cmpr, ig, olev, oc, usc, getmap, draw1D, draw2D):
    h5g = h5f['data']['grid_' + str(ig).zfill(10)]
    level = h5g.attrs['level']
    levok = (level == olev)
    if not levok:
        return levok, [], []

    ledge = h5g.attrs['left_edge']
    redge = h5g.attrs['right_edge']
    ngb = h5g.attrs['n_b']
    if oc is not None:
        inb, ind = pu.find_indices(ngb, oc, ledge, redge, draw1D, draw2D, False)
        if not any(inb):
            return False, [], []
        if not getmap:
            return levok, [[], inb, ledge / usc, redge / usc, olev, []], []
    off = h5g.attrs['off']
    n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
    ce = n_b + off
    dset = h5g[dset_name][:, :, :].swapaxes(0, 2)
    cmpr0, cmprb, h5c, cmprd, cmprl, cmprt, diff_struct = cmpr
    if cmpr0:
        dset = pu.execute_comparison(dset, h5c['data']['grid_' + str(ig).zfill(10)][cmprd][:, :, :].swapaxes(0, 2), cmprt)

    if oc is None:
        return pu.locate_extrema(dset, ledge, redge, ngb)

    b2d, b1d, extr = take_cuts_and_lines(dset, ind, draw1D, draw2D)

    return levok, [b2d, inb, ledge / usc, redge / usc, olev, b1d], extr


def take_cuts_and_lines(dset, ind, draw1D, draw2D):
    xy, xz, yz, d2min, d2max = [], [], [], [], []
    if draw2D[2]:
        xy = dset[:, :, ind[2]].swapaxes(0, 1)
        d2min.append(np.min(xy))
        d2max.append(np.max(xy))
    if draw2D[1]:
        xz = dset[:, ind[1], :].swapaxes(0, 1)
        d2min.append(np.min(xz))
        d2max.append(np.max(xz))
    if draw2D[0]:
        yz = dset[ind[0], :, :].swapaxes(0, 1)
        d2min.append(np.min(yz))
        d2max.append(np.max(yz))

    fx, fy, fz, d1min, d1max = [], [], [], [], []
    if draw1D[0]:
        fx = dset[:, ind[1], ind[2]]
        d1min.append(np.min(fx))
        d1max.append(np.max(fx))
    if draw1D[1]:
        fy = dset[ind[0], :, ind[2]]
        d1min.append(np.min(fy))
        d1max.append(np.max(fy))
    if draw1D[2]:
        fz = dset[ind[0], ind[1], :]
        d1min.append(np.min(fz))
        d1max.append(np.max(fz))

    d3min, d3max = np.min(dset), np.max(dset)
    if any(draw2D):
        d2max = max(d2max)
        d2min = min(d2min)
    if any(draw1D):
        d1max = max(d1max)
        d1min = min(d1min)

    return [yz, xz, xy], [fx, fy, fz], [d1min, d1max, d2min, d2max, d3min, d3max]


def collect_particles(h5f, drawh, center, player, uupd, usc, plotlevels, gridlist):
    if 'particle_types' not in list(h5f):
        return False, [], []
    print('Reading particles')
    px, py, pz, pm = np.array([]), np.array([]), np.array([]), np.array([])
    for gn in h5f['data']:
        if h5f['data'][gn].attrs['level'] in plotlevels and int(gn.split('grid_')[-1]) in gridlist:
            if str(player[1]) == '0' and str(player[2]) == '0' and str(player[3]) == '0':
                px = np.concatenate((px, h5f['data'][gn]['particles'][ps.particles_group]['position_x'][:]))
                py = np.concatenate((py, h5f['data'][gn]['particles'][ps.particles_group]['position_y'][:]))
                pz = np.concatenate((pz, h5f['data'][gn]['particles'][ps.particles_group]['position_z'][:]))
                if drawh:
                    pm = np.concatenate((pm, h5f['data'][gn]['particles'][ps.particles_group]['mass'][:]))
            else:
                apx = h5f['data'][gn]['particles'][ps.particles_group]['position_x'][:]
                apy = h5f['data'][gn]['particles'][ps.particles_group]['position_y'][:]
                apz = h5f['data'][gn]['particles'][ps.particles_group]['position_z'][:]
                maskx = np.abs(apx - center[0]) <= float(player[1])
                masky = np.abs(apy - center[1]) <= float(player[2])
                maskz = np.abs(apz - center[2]) <= float(player[3])
                if player[0]:
                    auxm = str(player[1]) == '0' or str(player[2]) == '0' or str(player[3]) == '0'
                    for i in range(np.size(maskx)):
                        maskx[i] = maskx[i] or masky[i] or maskz[i] or auxm
                else:
                    for i in range(np.size(maskx)):
                        maskx[i] = maskx[i] and masky[i] and maskz[i]
                px = np.concatenate((px, apx[maskx]))
                py = np.concatenate((py, apy[maskx]))
                pz = np.concatenate((pz, apz[maskx]))
                if drawh:
                    pm = np.concatenate((pm, h5f['data'][gn]['particles'][ps.particles_group]['mass'][maskx]))

    if px.size == 0:
        return False, [], []
    if uupd:
        return True, pu.list3_division([px, py, pz], usc), pm

    return True, [px, py, pz], pm
