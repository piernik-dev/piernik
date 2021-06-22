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
    refis, l2, h2, l3, h3 = [], [], [], [], []
    for iref in range(maxglev + 1):
        print('REFINEMENT ', iref)
        blks = []
        for ib in range(cgcount):
            block, extr = read_block(h5f, var, ib, iref, center, usc)
            blks.append(block)
            l2.append(extr[0])
            h2.append(extr[1])
            l3.append(extr[2])
            h3.append(extr[3])
        refis.append(blks)
    return refis, [min(l2), max(h2), min(l3), max(h3)]


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

    return [inb, [yz, xz, xy], ledge / usc, redge / usc], [d2min, d2max, d3min, d3max]


def collect_particles(h5f, nbins):
    print('Reading particles')
    px, py, pz, pm = np.array([]), np.array([]), np.array([]), np.array([])
    for gn in h5f['data']:
        px = np.concatenate((px, h5f['data'][gn]['particles']['stars']['position_x'][:]))
        py = np.concatenate((py, h5f['data'][gn]['particles']['stars']['position_y'][:]))
        pz = np.concatenate((pz, h5f['data'][gn]['particles']['stars']['position_z'][:]))
        if nbins > 1:
            pm = np.concatenate((pm, h5f['data'][gn]['particles']['stars']['mass'][:]))
    return [px, py, pz], pm
