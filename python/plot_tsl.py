#!/usr/bin/python

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse

remove_comments  = re.compile("(?!\#)", re.VERBOSE)

parser = argparse.ArgumentParser()
parser.add_argument("-f", nargs=1, default=None)
parser.add_argument("files", nargs='*')

args = parser.parse_args()
if len(args.files) < 1:
   parser.error("I need at least one tsl file")

data = []

for fn in args.files:
   f = open(fn,"rb")
   tab = [line.strip() for line in f.readlines()]
   f.close()
   header = np.array(tab[0][1:].split())
   if args.f == None:
      print ("There are following fields available in %s" % fn)
      print header
   else:
      field = args.f[0]
      fno = np.where(header == field)[0][0]

   tab = np.array([
      map(np.float64, line.split()) for line in filter(remove_comments.match, tab)
      ])
   data.append(tab)

fig = plt.figure()
ax = fig.add_subplot(111)
for i, fn in enumerate(data):
   ax.plot(fn[:, 1], fn[:, fno], label=args.files[i])

ax.legend()

plt.ylabel(field)
plt.xlabel(header[1])
plt.draw()
plt.show()
