#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, h5py, os, sys
import plot_compose as pc

if (len(sys.argv) < 3):
    print('PIERNIK VISUALIZATION FACILITY')
    print('Usage: ./pvf.py <file> <varname,[varname,...]> [options]')
    if len(sys.argv) < 2:
        exit()

cmap    = 'viridis'
plotdir = 'frames'
cu, cx, cy, cz = False, 0.0, 0.0, 0.0
zmin, zmax = 0.0, 0.0

def cli_params(argv):
    try:
        opts,args=getopt.getopt(argv,"c:ho:r:z:",["help","center=","colormap=","output=","zlim="])
    except getopt.GetoptError:
        print("Unidentified error.")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(" -h, \t\t--help \t\t\tprint this help \n\
 -c CX,CY,CZ, \t--center CX,CY,CZ \tplot cuts across given point coordinates CX, CY, CZ [default: computed domain center] \n\
 -o OUTPUT, \t--output OUTPUT \tdump plot files into OUTPUT directory [default: frames] \n\
 -r COLORMAP, \t--colormap COLORMAP \tuse COLORMAP palette [default: viridis] \n\
 -z ZMIN,ZMAX, \t--zlim ZMIN,ZMAX \tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]")
            sys.exit()

        elif opt in ("-c", "--center"):
            global cx, cy, cz, cu
            cx, cy, cz = arg.split(',')
            cu, cx, cy, cz = True, float(cx), float(cy), float(cz)

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

cli_params(sys.argv[3:])

pthfilen = sys.argv[1]
filen  = pthfilen.split('/')[-1]
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

h5f = h5py.File(pthfilen,'r')
if len(sys.argv) < 3:
    print("Available datafields: ", list(h5f['field_types'].keys()))
    exit(1)
if sys.argv[2] == "_all_":
    varlist = h5f['field_types'].keys()
else:
    varlist  = sys.argv[2].split(',')

print(varlist)

print("Reading file: %s" % pthfilen)
for var in varlist:
    #output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
    fnl = filen.split('/')[-1]
    output = plotdir+'/'+'_'.join(fnl.split('_')[:-1])+'_'+var+'_'+fnl.split('_')[-1].replace('.h5',".png")
    options = zmin, zmax, cmap, cu, cx, cy, cz
    pc.plotcompose(pthfilen, var, output, options)

h5f.close()
