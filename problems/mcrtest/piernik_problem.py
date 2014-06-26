#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('cairo')
from yt.mods import load as yt_load
from pylab import *

THRESHOLD = 1e-9
FIELD = "cr1"


def plot_diff(pf1, pf2, data1, data2, field):
    wd = pf1.domain_width
    n_d = pf1.domain_dimensions
    ext = np.array([pf1.domain_left_edge[1], pf1.domain_right_edge[1],
                    pf1.domain_left_edge[2], pf1.domain_right_edge[2]])
    ext *= pf1['pc']
    img1 = data1.to_frb(wd[0], (n_d[1] * 10, n_d[0] * 10),
                        center=np.array([0, 0, 0]), height=wd[1])
    img2 = data2.to_frb(wd[0], (n_d[1] * 10, n_d[0] * 10),
                        center=np.array([0, 0, 0]), height=wd[1])
    diff = (img2[field] - img1[field]) / (img1[field] + THRESHOLD)
    v = abs(diff).max()
    F = figure(1, (8, 6))
    imshow(diff, vmin=-v, vmax=v, extent=ext, cmap='RdBu')
    bar = colorbar()
    bar.ax.set_xlabel(
        r"$\frac{\rm{%s}^{\rm{new}} - \rm{%s}^{\rm{old}}}{\rm{%s}^{\rm{old}}}$" % (field, field, field))
    draw()
    xlabel('x [pc]')
    ylabel('y [pc]')
    savefig('diff.png')

    F = figure(1, (6, 6))
    clf()
    imshow(img2[field], extent=ext, cmap='algae')
    xlabel('x [pc]')
    ylabel('y [pc]')
    savefig('field.png')


if len(sys.argv) != 3:
    print("Wrong number of arguments!")
    sys.exit(-1)
PF1 = yt_load(sys.argv[1])
PF2 = yt_load(sys.argv[2])

DATA1 = PF1.h.slice(2, 0.0, fields=[FIELD])
DATA2 = PF2.h.slice(2, 0.0, fields=[FIELD])

if not PF1.h.field_list == PF2.h.field_list:
    print("Fields in files differ!")
    sys.exit(-1)

for field in PF1.h.field_list:
    if abs(DATA1[field] - DATA2[field]).max() >= THRESHOLD:
        print("Field %s differs" % field)
        plot_diff(PF1, PF2, DATA1, DATA2, field)
        sys.exit(-1)

figure(1, (8, 6))
draw()
savefig('diff.png')
savefig('field.png')
