#!/usr/bin/env python

from yt.mods import load
import sys
from matplotlib.pylab import imshow, savefig

for fn in sys.argv[1:]:
    fields = ['dend']

    pf = load(fn)
    c = 0.5 * (pf.domain_left_edge + pf.domain_right_edge)
    S = pf.domain_right_edge - pf.domain_left_edge
    n_d = pf.domain_dimensions

    slc = pf.h.slice(2, c[2], fields=fields)
    frb = slc.to_frb(S[0], (n_d[1], n_d[0]), height=S[1], center=c)
    imshow(frb['dend'])
    savefig('%s.png' % pf)
