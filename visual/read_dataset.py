#!/usr/bin/env python
import h5py as h5
import numpy as np


def collect_dataset(filen, dset_name):
    print('Reading', dset_name)
    h5f = h5.File(filen, 'r')
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

    h5f.close()
    return dset


def collect_particles(filen):
    print('Reading particles')
    px, py, pz = np.array([]), np.array([]), np.array([])
    h5f = h5.File(filen, 'r')
    for gn in h5f['data']:
        px = np.concatenate((px, h5f['data'][gn]['particles']['stars']['position_x'][:]))
        py = np.concatenate((py, h5f['data'][gn]['particles']['stars']['position_y'][:]))
        pz = np.concatenate((pz, h5f['data'][gn]['particles']['stars']['position_z'][:]))
    h5f.close()
    return px, py, pz


def convert_units(infile, toplot):
    au_cm = 1.49597870700e13
    pc_au = 206264.806248712
    pc_cm = pc_au * au_cm
    if infile == 'pc':
        cm = 1.0/pc_cm
        pc = 1.0
    elif infile == 'au':
        cm = 1.0/au_cm
        pc = pc_cm*cm
    elif infile == 'kpc':
        cm = 1.0/(1.0e3*pc_cm)
        pc = 0.001
    elif infile == 'm':
        cm = 1.0/1.0e2
        pc = pc_cm*cm
    else:
        return 1., False

    if toplot == 'cm':
        return cm, True
    elif toplot == 'metr':
        return 1.0e2*cm, True
    elif toplot == 'km':
        return 1.0e5*cm, True
    elif toplot == 'au':
        return au_cm*cm, True
    elif toplot == 'pc':
        return pc, True
    elif toplot == 'kpc':
        return 1.0e3*pc, True
    elif toplot == 'Mpc':
        return 1.0e6*pc, True
    elif toplot == 'lyr':
        return 9.4605e17*cm, True
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
