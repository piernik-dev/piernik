#!/usr/bin/env python
import h5py as h5
import numpy as np
import plot_utils as pu


def reconstruct_uniform(h5f, var, cu, center, nd, smin, smax):
    dset = collect_dataset(h5f, var)

    if cu:
        inb, ind = pu.find_indices(nd, center, smin, smax, True)
        print('Ordered plot center', center[0], center[1], center[2], ' gives following uniform grid indices:', ind[0], ind[1], ind[2])
    else:
        ind = int(nd[0] / 2), int(nd[1] / 2), int(nd[2] / 2)

    xy = dset[:, :, ind[2]]
    xz = dset[:, ind[1], :].swapaxes(0, 1)
    yz = dset[ind[0], :, :].swapaxes(0, 1)

    d3min, d3max = np.min(dset), np.max(dset)
    d2max = max(np.max(xz), np.max(xy), np.max(yz))
    d2min = min(np.min(xz), np.min(xy), np.min(yz))
    return xy, xz, yz, [d2min, d2max, d3min, d3max]


def collect_dataset(h5f, dset_name):
    print('Reading', dset_name)
    attrs = h5f['domains']['base'].attrs
    nxd, nyd, nzd = attrs['n_d'][0:3]
    dset = np.zeros((nxd, nyd, nzd))

    print('Reconstructing domain from cg parts')
    grid = h5f['grid_dimensions']
    for ig in range(grid.shape[0]):
        h5g = h5f['data']['grid_' + str(ig).zfill(10)]
        if h5g.attrs['level'] == 0:
            off = h5g.attrs['off']
            ngb = h5g.attrs['n_b']
            n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
            ce = n_b + off
            dset[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] = h5g[dset_name][:, :, :].swapaxes(0, 2)

    return dset


def collect_gridlevels(h5f, var, maxglev, cgcount, center, usc):
    refis = []
    for iref in range(maxglev + 1):
        print('REFINEMENT ', iref)
        blks = []
        for ib in range(cgcount):
            block = read_block(h5f, var, ib, iref, center, usc)
            blks.append(block)
        refis.append(blks)
    return refis


def read_block(h5f, dset_name, ig, olev, oc, usc):
    h5g = h5f['data']['grid_' + str(ig).zfill(10)]
    level = h5g.attrs['level']
    ledge = h5g.attrs['left_edge']
    redge = h5g.attrs['right_edge']
    ngb = h5g.attrs['n_b']
    levok = (level == olev)
    inb, ind = pu.find_indices(ngb, oc, ledge, redge, False)

    clen = h5g.attrs['dl']
    off = h5g.attrs['off']
    n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
    ce = n_b + off
    dset = h5g[dset_name][:, :, :].swapaxes(0, 2)
    xy = dset[:, :, ind[2]]
    xz = dset[:, ind[1], :].swapaxes(0, 1)
    yz = dset[ind[0], :, :].swapaxes(0, 1)
    d3min, d3max = np.min(dset), np.max(dset)
    d2max = max(np.max(xz), np.max(xy), np.max(yz))
    d2min = min(np.min(xz), np.min(xy), np.min(yz))

    return inb, xy, xz, yz, ledge / usc, redge / usc


def collect_particles(h5f, nbins):
    print('Reading particles')
    px, py, pz, pm = np.array([]), np.array([]), np.array([]), np.array([])
    for gn in h5f['data']:
        px = np.concatenate((px, h5f['data'][gn]['particles']['stars']['position_x'][:]))
        py = np.concatenate((py, h5f['data'][gn]['particles']['stars']['position_y'][:]))
        pz = np.concatenate((pz, h5f['data'][gn]['particles']['stars']['position_z'][:]))
        if nbins > 1:
            pm = np.concatenate((pm, h5f['data'][gn]['particles']['stars']['mass'][:]))
    return px, py, pz, pm


def convert_units(infile, toplot):
    au_cm = 1.49597870700e13
    pc_au = 206264.806248712
    pc_cm = pc_au * au_cm
    if infile == 'pc':
        cm = 1.0 / pc_cm
        pc = 1.0
    elif infile == 'au':
        cm = 1.0 / au_cm
        pc = pc_cm * cm
    elif infile == 'kpc':
        cm = 1.0 / (1.0e3 * pc_cm)
        pc = 0.001
    elif infile == 'm':
        cm = 1.0 / 1.0e2
        pc = pc_cm * cm
    else:
        return 1., False

    if toplot == 'cm':
        return cm, True
    elif toplot == 'metr':
        return 1.0e2 * cm, True
    elif toplot == 'km':
        return 1.0e5 * cm, True
    elif toplot == 'au':
        return au_cm * cm, True
    elif toplot == 'pc':
        return pc, True
    elif toplot == 'kpc':
        return 1.0e3 * pc, True
    elif toplot == 'Mpc':
        return 1.0e6 * pc, True
    elif toplot == 'lyr':
        return 9.4605e17 * cm, True
    else:
        return 1., False


def change_units(fromfile, toplot):
    infile = fromfile.decode('utf-8')
    if infile == toplot or toplot == '':
        return 1., fromfile, False
    if toplot == 'k':
        return 1.e3, b"".join([b'k', fromfile]), True
    if toplot == 'M':
        return 1.e6, b"".join([b'M', fromfile]), True
    if toplot == 'G':
        return 1.e9, b"".join([b'G', fromfile]), True
    if toplot == 'm':
        return 1.e-3, b"".join([b'm', fromfile]), True
    if toplot == 'mu':
        return 1.e-6, b"".join([b'mu', fromfile]), True
    conv, chan = convert_units(infile, toplot)
    if chan:
        return conv, bytes(toplot, 'utf-8'), chan
    else:
        return conv, fromfile, chan
