#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt
import h5py
import os
import sys
import plot_compose as pc

cmap = 'viridis'
plotdir = 'frames'
sctype = 'linear'
cu, cx, cy, cz = False, 0.0, 0.0, 0.0
zmin, zmax = 0.0, 0.0

print('PIERNIK VISUALIZATION FACILITY')

def print_usage():
    print('Usage: ./pvf.py <file> <varname,[varname,...] | _all_> [options]')
    print('')
    print('./pvf.py <file> - print available datafields from <file>')
    print('./pvf.py <file> _all_ - plot all available datafields from <file>')
    print('./pfv.py <file> <varname,[varname,...] [options] - plot specified datafields from <file>')
    print('')
    print('Options:')
    print(' -h, \t\t--help \t\t\tprint this help')
    print(' -c CX,CY,CZ, \t--center CX,CY,CZ \tplot cuts across given point coordinates CX, CY, CZ [default: computed domain center]')
    print(' -l SCALETYPE, \t--scale SCALETYPE \tdump use SCALETYPE scale type for displaying data (possible values: 0 | linear, 1 | symlin, 2 | log, 3 | symlog) [default: linear]')
    print(' -o OUTPUT, \t--output OUTPUT \tdump plot files into OUTPUT directory [default: frames]')
    print(' -r COLORMAP, \t--colormap COLORMAP \tuse COLORMAP palette [default: viridis]')
    print(' -z ZMIN,ZMAX, \t--zlim ZMIN,ZMAX \tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]')

def cli_params(argv):
    try:
        opts, args = getopt.getopt(argv, "c:hl:o:r:z:", ["help", "center=", "colormap=", "output=", "scale=", "zlim="])
    except getopt.GetoptError:
        print("Unrecognized options: %s \n" % argv)
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print_usage()
            sys.exit()

        elif opt in ("-c", "--center"):
            global cx, cy, cz, cu
            cx, cy, cz = arg.split(',')
            cu, cx, cy, cz = True, float(cx), float(cy), float(cz)

        elif opt in ("-l", "--scale"):
            global sctype
            sctype = str(arg)

        elif opt in ("-o", "--output"):
            global plotdir
            plotdir = str(arg)
            print('PLOTDIR: ', plotdir)

        elif opt in ("-r", "--colormap"):
            global cmap
            cmap = str(arg)

        elif opt in ("-z", "--zlim"):
            global zmin, zmax
            zmin, zmax = arg.split(',')
            zmin = float(zmin)
            zmax = float(zmax)
            print("zmin, zmax = ", zmin, zmax)

def check_file(pfile):
    if not os.path.exists(pfile):
        print('The file %s does not exist!' % pfile)
        exit()

def list_file_fields(pfile):
    h5f = h5py.File(pfile, 'r')
    print('Available datafields in the file %s: \n' % pfile, list(h5f['field_types'].keys()))

pthfilen = sys.argv[1]
if (len(sys.argv) < 3):
    print_usage()
    if (len(sys.argv) == 2 and sys.argv[-1][0] != '-'):
        print('')
        check_file(pthfilen)
        list_file_fields(pthfilen)
        exit()

cli_params(sys.argv[3:])

check_file(pthfilen)
h5f = h5py.File(pthfilen, 'r')

filen = pthfilen.split('/')[-1]
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

if sys.argv[2] == "_all_":
    varlist = h5f['field_types'].keys()
else:
    varlist = sys.argv[2].split(',')

print(varlist)

print("Reading file: %s" % pthfilen)
for var in varlist:
    # output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
    fnl = filen.split('/')[-1]
    output = plotdir + '/' + '_'.join(fnl.split('_')[:-1]) + '_' + var + '_' + fnl.split('_')[-1].replace('.h5', ".png")
    options = zmin, zmax, cmap, sctype, cu, cx, cy, cz
    pc.plotcompose(pthfilen, var, output, options)

h5f.close()
