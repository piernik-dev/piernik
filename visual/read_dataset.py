#!/usr/bin/env python


def collect_dataset(filen, dset_name):
    import h5py as h5
    import numpy as np
    print('Prepearing to read ', dset_name)
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
