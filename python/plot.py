#!/usr/bin/python
import sys
import tables as h5
from optparse import OptionParser
from pylab import *

def read_hdf5(options,filename):
   h5f = h5.openFile(filename,'r')
   dsets = [node.name for node in h5f.walkNodes("/","Array")]
   attrs = h5f.root._v_attrs

   if not options.dset in dsets:
      print '\033[91m' + "unknown dataset: " + '\033[0m' + options.dset
      print 'Possible choices are: '
      print dsets
      exit(-1)

   arr = h5f.getNode("/",name=options.dset)[:,:,:]

   h5f.close()
   return arr, attrs

usage = "usage: %prog [options] FILE"
try:
   parser = OptionParser(usage=usage, epilog="if you want some additional features, bug me.")
except TypeError:
   parser = OptionParser(usage=usage)

parser.add_option("-d", "--dataset", dest="dset", default="",
      help="datasets, use comma-separated list")
parser.add_option("-p", "--plane", dest="plane", default="xy",
      help="plane, possible values {xy,yz,xz}")
parser.add_option("-c", "--cell", dest="cell", default=0, type="int",
      help="plane, possible values {xy,yz,xz}")
parser.add_option("--log", dest="dolog", default=False, action="store_true",
      help="use logarithmic scal")
(options, args) = parser.parse_args(sys.argv[1:])

if not options.plane in ["xy","yz","xz"]:
   print '\033[91m' + "wrong plane: " + '\033[0m' + options.plane
   exit(-1)

filename = args[0]
tab, attrs = read_hdf5(options,filename)
shape = list(tab.shape)

if options.plane == "xy":
   p = 0
   extent = [attrs.xmin[0], attrs.xmax[0], attrs.ymin[0], attrs.ymax[0]]
   s = s_[:,:,options.cell]
elif options.plane == "xz":
   p = 1
   extent = [attrs.xmin[0], attrs.xmax[0], attrs.zmin[0], attrs.zmax[0]]
   s = s_[:,options.cell,:]
else:
   p = 2
   extent = [attrs.ymin[0], attrs.ymax[0], attrs.zmin[0], attrs.zmax[0]]
   s = s_[:,:,options.cell]

if options.cell >= shape[p] or options.cell < 0:
   print("You gotta be kiddin' me... %d is _not_ in [0,%d]" % (options.cell, shape[p]-1))

arr = tab[s]
if options.dolog:
   arr = np.log10(arr)

title_str = "T = %f" % attrs.time[0]

imshow(arr,extent=extent)
title(title_str)
colorbar()
draw()
show()
