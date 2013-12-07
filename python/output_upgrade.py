#!/usr/bin/python
import h5py as h5
import numpy as np
import sys


def upgrade_200_201(fname):
    with h5.File(fname, 'r+') as h5f:
        ver = h5f.attrs['piernik']
        if (ver < 2.01) and (ver >= 2.00):
            print("\t|- Applying fix for 2.00 -> 2.01")
            for key, group in h5f['field_types'].iteritems():
                if 'field_to_cgs' in group.attrs.keys():
                    group.attrs['field_to_cgs'] = \
                        1.0 / group.attrs['field_to_cgs']
            h5f.attrs['piernik'] = np.array([2.01])

for fname in sys.argv[1:]:
    try:
        with h5.File(fname, 'r') as h5f:
            try:
                ver = h5f.attrs['piernik'][0]
            except AttributeError:
                print("%s is not PIERNIK output" % fname)
                continue
    except IOError:
        print("Cannot open %s. It's either not hdf5 or doesn't exist" % fname)
        continue
    print("Scouring file %s" % fname)
    upgrade_200_201(fname)
    print("\t|___ DONE")
