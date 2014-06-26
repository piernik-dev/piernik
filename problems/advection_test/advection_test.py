#!/usr/bin/env python

from yt.config import ytcfg
ytcfg["yt", "loglevel"] = "50"
from yt.mods import load
import numpy as np


def calculate_norm(fn):
    pf = load(fn)

    data = pf.h.all_data()
    diff = data['inid'] - data['denn']
    norm = np.sqrt((diff ** 2).sum() / data['inid'].sum())
    print("Calculating difference between numerical and analytical solution")
    print("L2 error norm = %12.6f, min and max error = %15.6f %15.6f" %
          (norm, diff.min(), diff.max()))
    del data, diff, pf


def diff_files(fn1, fn2, field):
    pf1 = load(fn1)
    data1 = pf1.h.all_data()[field]
    pf2 = load(fn2)
    data2 = pf2.h.all_data()[field]
    diff = data1 - data2
    norm = np.sqrt((diff ** 2).sum() / data1.sum())
    print("Calculating difference in %s between files %s and %s" %
          (field, pf1, pf2))
    print("L2 error norm = %12.6f, min and max error = %15.6f %15.6f" %
          (norm, diff.min(), diff.max()))
    del data1, data2, diff, pf1, pf2

for i in range(1, 3):
    fn1 = "moving_pulse_ts1_%4.4i.h5" % i
    fn2 = "moving_pulse_ts2_%4.4i.h5" % i
    calculate_norm(fn1)
    calculate_norm(fn2)
    diff_files(fn1, fn2, 'denn')
