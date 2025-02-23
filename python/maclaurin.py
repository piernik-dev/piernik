#!/usr/bin/env python

# This file requires use of Piernik I/O format v 1.x

import os


def Maclaurin_test(file):

    missing = []
    try:
        import h5py as h5
    except ImportError:
        missing.append("h5py")
    try:
        import numpy as np
    except ImportError:
        missing.append("NumPy")
    try:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import NullFormatter
    except ImportError:
        missing.append("matplotlib")

    if (len(missing) > 0):
        print("You must install the package(s) ", missing)
        return

    # matplotlib.use('cairo')      # choose output format
    # This may sometimes help with font issues
    # from matplotlib import rc
    # rc('text',usetex=True)

    if (not os.path.isfile(file)):
        print("Cannot find ", file)
        return

    try:
        h5f = h5.File(file, "r")
    except IOError:
        print("Cannot open '" + file + "' as HDF5")
        return

    xmin, xmax = h5f['domains/base'].attrs['x-edge_position']
    ymin, ymax = h5f['domains/base'].attrs['y-edge_position']
    zmin, zmax = h5f['domains/base'].attrs['z-edge_position']

    fpiG = h5f.attrs['fpiG']
    a = h5f.attrs['a1']
    x0 = h5f.attrs['x0'][0]
    y0 = h5f.attrs['y0'][0]
    z0 = h5f.attrs['z0'][0]

    n = 0
    for dname in h5f['data'].keys():
        lev_max = max(0, h5f['data'][dname].attrs['level'])
        n += int(np.prod(h5f['data'][dname].attrs['n_b']))

    nz, ny, nx = h5f['domains/base'].attrs['n_d'] * 2**lev_max
    dx, dy, dz = (xmax - xmin) / nx, (ymax - ymin) / ny, (zmax - zmin) / nz

    r = np.zeros(n)
    phi = np.zeros(n)
    phi0 = np.zeros(n)

    # *********************************************************************

    ind = 0
    for dname in h5f['/data'].keys():
        try:
            soln = h5f['data'][dname]['gpot'][:, :, :]
            phi0_3D = h5f['data'][dname]['apot'][:, :, :]
        except KeyError:
            print("Cannot find apot or gpot arrays")
            return

        off = h5f['data'][dname].attrs['off']
        lev = h5f['data'][dname].attrs['level'][0]

        sc_fac = 2. ** (lev_max - lev)
        for i in range(0, int(h5f['data'][dname].attrs['n_b'][0])):
            x2 = (xmin + (off[0] + i + 0.5) * sc_fac * dx - x0)**2

            for j in range(0, int(h5f['data'][dname].attrs['n_b'][1])):
                y2 = (ymin + (off[1] + j + 0.5) * sc_fac * dy - y0)**2

                for k in range(0, int(h5f['data'][dname].attrs['n_b'][2])):
                    z2 = (zmin + (off[2] + k + 0.5) * sc_fac * dz - z0)**2

                    r[ind] = np.sqrt(x2 + y2 + z2)
                    phi[ind] = soln[k, j, i]
                    phi0[ind] = phi0_3D[k, j, i]
                    ind += 1

    h5f.close()

    new_r = np.unique(r)
    mean = np.zeros_like(new_r)
    std = np.zeros_like(new_r)
    phi_0 = np.zeros_like(new_r)

    for i in range(0, new_r.shape[0] - 1):
        ind = np.where(r == new_r[i])
        temp = phi[ind]
        temp2 = phi0[ind]
        phi_0[i] = temp2.mean()
        mean[i] = temp.mean()
        std[i] = temp.std()
    # *********************************************************************
    # plotting ---------------------------------
    # definitions for the axes
    left, width = 0.15, 0.80

    rect_hi = [left, 0.3, width, 0.65]
    rect_lo = [left, 0.05, width, 0.2]

    fig = plt.figure(1, figsize=(8, 8))

    axhi = fig.add_axes(rect_hi)
    axlo = fig.add_axes(rect_lo)

    GM = fpiG / 3.

    axhi.set_title("Maclaurin spheroid $e=0, a_1=1, \\varrho_0=1$")
    axhi.plot(new_r[1:-1], mean[1:-1] / GM,
              'g.', new_r[1:-1], phi_0[1:-1] / GM, 'b')
    axhi.xaxis.set_major_formatter(NullFormatter())
    axhi.set_ylabel('Gravitational potential / GM')
    axhi.legend(('Numerical solution - $\\varphi$',
                 'Analytical solution - $\\varphi_0$'), loc='lower right')

    # Typically, the lower right position is best for the legend because at
    # lmax=16 the multipole solver underestimates the image mass.
    # For lmax near 192 the image mass almost matches the mass of the sphere.
    # For lmax much greates than 192 it may help a bit to set loc to 'upper
    # right'
    axlo.plot(new_r[1:-1], (phi_0[1:-1] - mean[1:-1]) / GM, 'g', new_r[1:-1],
              (phi_0[1:-1] - mean[1:-1] - std[1:-1]) / GM, 'r:', new_r[1:-1],
              (phi_0[1:-1] - mean[1:-1] + std[1:-1]) / GM, 'r:')
    axlo.legend(('avg. difference', '+/- deviation'), loc='lower right')
    axlo.set_xlabel('Radius')
    axlo.set_ylabel('($\\varphi_0 - \\varphi$) / GM')

    plt.savefig('maclaurin.png', facecolor='white')
    plt.close(fig)


if __name__ == "__main__":
    import sys
    if (len(sys.argv) <= 1):
        print("Usage : ", sys.argv[0], " maclaurin_test_hdf_file")
    else:
        Maclaurin_test(sys.argv[1])
        if (len(sys.argv) > 2):
            print("Ignored arguments: ", sys.argv[2:len(sys.argv)])
