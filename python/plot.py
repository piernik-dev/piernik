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

def get_frame(filename,options):
   tab, attrs = read_hdf5(options,filename)
   shape = list(tab.shape)
   i1 = options.i1
   i2 = options.i2

   if options.plane == "xy":
      p = 0
      extent = [attrs.xmin[0], attrs.xmax[0], attrs.ymin[0], attrs.ymax[0]]
      s = s_[:,:,i1]
   elif options.plane == "xz":
      p = 1
      extent = [attrs.xmin[0], attrs.xmax[0], attrs.zmin[0], attrs.zmax[0]]
      s = s_[:,i1,:]
   elif options.plane == "yz":
      p = 2
      extent = [attrs.ymin[0], attrs.ymax[0], attrs.zmin[0], attrs.zmax[0]]
      s = s_[:,:,i1]
   elif options.plane == "x":
      p = 0
      extent = [attrs.xmin[0], attrs.xmax[0]]
      s = s_[i2,i1,:]
   elif options.plane == "y":
      p = 1
      extent = [attrs.ymin[0], attrs.ymax[0]]
      s = s_[i1,:,i2]
   elif options.plane == "z":
      p = 2
      extent = [attrs.zmin[0], attrs.zmax[0]]
      s = s_[:,i2,i1]

   if i1 >= shape[p] or i1 < 0:
      print("You gotta be kiddin' me... %d is _not_ in [0,%d]" % (i1, shape[p]-1))

   arr = tab[s]
   if options.dolog:
      arr = np.log10(arr)

   title_str = "T = %f" % attrs.time[0]

   return arr, title_str, extent

usage = "usage: %prog [options] FILE"
try:
   parser = OptionParser(usage=usage, epilog="if you want some additional features, bug me.")
except TypeError:
   parser = OptionParser(usage=usage)

parser.add_option("-d", "--dataset", dest="dset", default="",
      help="datasets, use comma-separated list")
parser.add_option("-p", "--plane", dest="plane", default="xy",
      help="plane, possible values {xy,yz,xz}")
parser.add_option("-c", "--cell", "--i1", dest="i1", default=0, type="int",
      help="")
parser.add_option("--i2", dest="i2", default=0, type="int",
      help="")
parser.add_option("--log", dest="dolog", default=False, action="store_true",
      help="use logarithmic scal")
(options, args) = parser.parse_args(sys.argv[1:])

if not options.plane in ["xy","yz","xz","x","y","z"]:
   print '\033[91m' + "wrong plane: " + '\033[0m' + options.plane
   exit(-1)

if len(args) > 1:
   ion()

arr, title_str, extent = get_frame(args[0], options)
if len(extent) == 2:
   plt = plot(linspace(extent[0],extent[1],num=len(arr)), arr)
else:
   plt = imshow(arr,extent=extent,interpolation='nearest')
   col=colorbar()
title(title_str)
draw()
show()

if len(args) > 1:
   for arg in args[1:]:
      arr, title_str, extent = get_frame(arg,options)
      if len(extent) == 2:
         plt = plot(linspace(extent[0],extent[1],num=len(arr)), arr)
      else:
         plt = imshow(arr,extent=extent,interpolation='nearest')
         col.update_bruteforce(plt)
      title(title_str)
      draw()
   print "Press ENTER to exit"
   e = sys.stdin.readline()  # wait for ENTER due to ion()
