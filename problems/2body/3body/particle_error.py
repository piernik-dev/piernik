#!/usr/bin/env python3
import sys
import h5py
import numpy as np
from numpy.linalg import norm
from math import sqrt

pkeys = ['id', 'mass', 'energy', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z', 'acceleration_x', 'acceleration_y', 'acceleration_z']


def pr_err(v):
    print("3-body failed")
    print("%f" % v)
    exit(0)


if (len(sys.argv) != 2):
    sys.stderr.write("Error: only two arguments accepted.\nUsage: " + sys.argv[0] + " piernik_hdf_file\n")
    exit(1)

try:
    h5f = h5py.File(sys.argv[1], "r")
except:
    pr_err(.3)

try:
    n = h5f['grid_particle_count'][0]
except:
    pr_err(.2)

# This is broken when file came from a parallel run
# if n != 3:
#     pr_err(.1)

m = [0., 0., 0.]
L = [0., 0., 0.]
for ig in range(h5f['grid_dimensions'].shape[0]):
    h5g = h5f['data']['grid_' + str(ig).zfill(10)]
    try:
        h5gps = h5g["particles"]["stars"]
        for i in range(h5gps.attrs['n_part'][0]):
            x = [h5gps['position_x'][i], h5gps['position_y'][i], h5gps['position_z'][i]]
            v = [h5gps['velocity_x'][i], h5gps['velocity_y'][i], h5gps['velocity_z'][i]]
            if h5gps['id'][i] == 3:  # hardcoded

                # O-C or period error for the central body
                dt = - np.dot(x, v) / np.dot(v, v)

                # projected closest approach – not a good estimate of error
                # if we really want it, we need to estimate it from two closest points during approach to the origin (in-code measurement)
                # b = sqrt(cdot(x, x) - abs(cdot(x, v) * dt))

                # other nice error estimates:
                # * angle between velocity at closest approach and initial velocity (requires in-code measurement)

            if h5gps['id'][i] == 1:  # hardcoded
                m1 = np.multiply(x, v)
                L1 = np.cross(x, v)

            # total momentum w.r.t. momentum of particle 1 or 2
            m = np.add(m, v)

            # total angular momentum w.r.t. angular momentum of particle 1 or 2
            L += np.cross(x, v)

            # todo: change of total energy (requires initial values)
    except:
        pr_err(.4)

print("Period error,momentum error,angular momentum error")
print("%g , %g , %g" % ((0. if dt == 0. else abs(dt / h5f.attrs["time"][0])), norm(m) / norm(m1), norm(L) / norm(L1)))
