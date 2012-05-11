#/usr/bin/env python

import sys
import h5py as h5
import numpy as np

def die(h5f, msg):
    h5f.close()
    print("%s" % msg)
    sys.exit()

def grow_dsets(h5f, datasets):
    for dset in datasets:
        fluid = h5f[dset][:].T.copy()
        n, nx, ny, nz = fluid.shape
        del h5f[dset]
        h5f[dset] = np.resize(fluid, (n+1, nx, ny, nz)).T
        h5f[dset][:,:,:,-1] = 0.0

for f in sys.argv[1:]:
    h5f = h5.File(f)
    try:
       print("Piernik Version %f"%(h5f.attrs['piernik']))
    except:
       die("%s is not piernik restart!!!"%(f))

    if h5f.attrs['piernik'][0] >= 2.0:
        for grid in h5f['data'].keys():
            grow_dsets(h5f['data'][grid], ['fluid', 'u_0'])
    else:
        grow_dsets(h5f, ['fluid', 'u_0'])

    h5f.close()
