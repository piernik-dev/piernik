#!/usr/bin/env python

import tables as h5
import numpy as np
from pylab import *
from optparse import OptionParser
import time,sys

def read_data_hdf5(name):
   h5f = h5.openFile(name)
   attrs = (h5f.root._v_attrs.xmin[0], h5f.root._v_attrs.xmax[0], h5f.root._v_attrs.time[0])
   denn = h5f.root.denn[0,:,:]
   vlxn = h5f.root.vlxn[0,:,:]
   vlyn = h5f.root.vlyn[0,:,:]
   dend = h5f.root.dend[0,:,:]
   vlxd = h5f.root.vlxd[0,:,:]
   vlyd = h5f.root.vlyd[0,:,:]
   h5f.close()
   return denn, dend, vlxn, vlxd, vlyn, vlyd, attrs

def sum_y(tab,n):
   return np.sum(tab.T,1)/n

usage = "usage: %prog FILE"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()
if len(args) < 1:
      parser.error("I need at least one file")

plt.ion()

fig = plt.figure(1, figsize=(10,12))

denn, dend, vlxn, vlxd, vlyn, vlyd, attrs = read_data_hdf5(args[0])
ny,nx = denn.shape
xmin = attrs[0]
xmax = attrs[1]
time = attrs[2]
x = np.linspace(xmin,xmax,nx)

subplot(211)
axis([xmin,xmax,0.01,1.0])
loglog(x, sum_y(denn,ny), label='denn_0')
loglog(x, sum_y(dend,ny), label='dend_0')
loglog(x, sum_y(vlyn,ny), label='vlyn_0')
loglog(x, sum_y(vlyd,ny), label='vlyd_0')
sm1, = loglog(x,sum_y(denn,ny))
sm2, = loglog(x,sum_y(dend,ny))
sm3, = loglog(x,sum_y(vlyn,ny))
sm4, = loglog(x,sum_y(vlyd,ny))
legend(loc=3)

subplot(212)
axis([xmin,xmax,-0.001,0.001])
p1, = plot(x, sum_y(vlxd,ny), label='vlxd')
p2, = plot(x, sum_y(vlxn,ny), label='vlxn')
legend()
draw()

if len(args) > 1:
   for name in args:
      denn, dend, vlxn, vlxd, vlyn, vlyd, attrs = read_data_hdf5(name)
      xmin = attrs[0]
      xmax = attrs[1]
      time = attrs[2]
      subplot(211)
      sm1.set_ydata(sum_y(denn,ny))
      sm1.set_label("denn")
      sm2.set_ydata(sum_y(dend,ny))
      sm2.set_label("dend")
      sm3.set_ydata(sum_y(vlyn,ny))
      sm3.set_label("vlyn")
      sm4.set_ydata(sum_y(vlyd,ny))
      sm4.set_label("vlyd")
      sm1.axes.set_title("T = %f" % time)
      legend(loc=3)
      subplot(212)
      p1.set_ydata(sum_y(vlxd,ny))
      p2.set_ydata(sum_y(vlxn,ny))
      draw()

e = sys.stdin.readline()  # czekaj na enter
