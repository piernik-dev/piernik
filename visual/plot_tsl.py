#!/usr/bin/env python

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse

remove_comments = re.compile("(?!\#)", re.VERBOSE)

parser = argparse.ArgumentParser()
parser.add_argument("-f", nargs=1, default=None, help='variable to plot on vertical axis')
parser.add_argument("-x", nargs=1, default=['time'], help='variable to plot on horizontal axis')
parser.add_argument("--xlim", nargs=2, default=None, help='scale span for horizontal axis')
parser.add_argument("--ylim", nargs=2, default=None, help='scale span for vertical axis')
parser.add_argument("-s", "--save", nargs=1, default=None, help="save plot to a file")
parser.add_argument("files", nargs='*')

args = parser.parse_args()
if len(args.files) < 1:
    parser.error("I need at least one tsl file")

data = []
nenough = True


def print_header(fn, header):
    print("There are following fields available in %s" % fn)
    for ii, he in enumerate(header):
        print(ii+1, he)


for fn in args.files:
    f = open(fn, "r")
    tab = [line.strip() for line in f.readlines()]
    f.close()
    header = np.array(tab[0][1:].split())

    if args.f is None:
        nenough = True
    elif args.f[0] not in header:
        print("Field %s not met in %s!" % (args.f[0], fn))
    elif args.x[0] not in header:
        print("Field %s not met in %s!" % (args.x[0], fn))
    else:
        nenough = False
    if nenough:
        print_header(fn, header)
        exit(0)

    field = args.f[0]
    fno = np.where(header == field)[0][0]
    xtime = args.x[0]
    xno = np.where(header == xtime)[0][0]

    tab_work = tab   # problems with map in python3+
    tab = np.array([[float(item) for item in line.split()]
                    for line in filter(remove_comments.match, tab_work)])
    data.append(tab)

fig = plt.figure()
ax = fig.add_subplot(111)
if args.xlim is not None:
    plt.xlim(float(args.xlim[0]), float(args.xlim[1]))
if args.ylim is not None:
    plt.ylim(float(args.ylim[0]), float(args.ylim[1]))

for i, fn in enumerate(data):
    ax.plot(fn[:, xno], fn[:, fno], label=args.files[i])

ax.legend()

plt.ylabel(field)
plt.xlabel(xtime)
plt.draw()
if args.save is None:
    plt.show()
else:
    plt.savefig(args.save[0], facecolor='white')
    print(args.save[0], "written to disk")
plt.clf()
