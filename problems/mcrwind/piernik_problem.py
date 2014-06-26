#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('cairo')
from yt.mods import load as yt_load
from pylab import *

THRESHOLD = 1e-9
FIELD = "cr1"

def _myplot(diff, fname, ext, clbl):
    v = abs(diff).max()
    figure(1, (6, 8))
    imshow(diff, vmin=-v, vmax=v, extent=ext, cmap='RdBu')
    bar = colorbar()
    bar.ax.set_xlabel(clbl)
    draw()
    xlabel('y [pc]')
    ylabel('z [pc]')
    savefig(fname)
    clf()

def plot_diff(pf1, pf2, data1, data2, field):
    wd = pf1.domain_width
    n_d = pf1.domain_dimensions
    ext = np.array([pf1.domain_left_edge[1], pf1.domain_right_edge[1],
                    pf1.domain_left_edge[2], pf1.domain_right_edge[2]])
    ext *= pf1['pc']
    img1 = data1.to_frb(wd[1], (n_d[2] * 10, n_d[1] * 10),
                        center=np.array([0, 0, 0]), height=wd[2])
    img2 = data2.to_frb(wd[1], (n_d[2] * 10, n_d[1] * 10),
                        center=np.array([0, 0, 0]), height=wd[2])
    diff = (img2[field] - img1[field])
    clbl = \
        r"$\rm{%s}^{\rm{new}} - \rm{%s}^{\rm{old}}$" % (field, field)
    _myplot(diff, 'diff_bare.png', ext, clbl)

    clbl = \
        r"$\frac{\rm{%s}^{\rm{new}} - \rm{%s}^{\rm{old}}}{\rm{%s}^{\rm{old}}}$" % (field, field, field)
    _myplot(diff / (img1[field] + THRESHOLD), 'diff.png', ext, clbl)


if len(sys.argv) != 3:
    print("Wrong number of arguments!")
    sys.exit(-1)
PF1 = yt_load(sys.argv[1])
PF2 = yt_load(sys.argv[2])

axis = np.where(PF1.h.grids[0].ActiveDimensions == 1)[0][0]
DATA1 = PF1.h.slice(axis, 0.0, fields=[FIELD])
DATA2 = PF2.h.slice(axis, 0.0, fields=[FIELD])

if not PF1.h.field_list == PF2.h.field_list:
    print("Fields in files differ!")
    sys.exit(-1)

for field in PF1.h.field_list:
    if abs(DATA1[field] - DATA2[field]).max() >= THRESHOLD:
        print("Field %s differs" % field)
        plot_diff(PF1, PF2, DATA1, DATA2, field)
        sys.exit(-1)

figure(1, (8,6))
draw()
savefig('diff.png')
savefig('diff_bare.png')
