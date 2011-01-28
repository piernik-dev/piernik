#!/usr/bin/python

import sys
import numpy as np
import pylab as P
from optparse import OptionParser

usage = "usage: %prog FILE COLUMN_NUMBER"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()
if len(args) < 1:
   parser.error("I need at least one file and column number")

f = open(args[0],"rb")
tab = []
if len(args) < 2:
   line = f.readline().strip().split()
else:
   for line in f.readlines():
      tab.append(line.strip().split())

f.close()
if (len(tab) < 1):
   print ("Data available in %s file." % args[0])
   for i in range(0,len(line)):
      print (" %2i - %s " % (i,line[i]))
   exit()

names = tab.pop(0)
tab.remove(["#"])

data=np.array(tab)

plt = P.plot(data[:,1], data[:,args[1]])
P.ylabel(names[int(args[1])])
P.xlabel(names[1])
P.show()
