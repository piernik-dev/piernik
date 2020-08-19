#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import h5py
import plot_compose as pc

if (len(sys.argv) < 3):
    print('PIERNIK VISUALIZATION FACILITY')
    print('Usage: ./pvf.py <file> <varname [varname ...]>')
    if len(sys.argv) < 2:
        exit()

pthfilen = sys.argv[1]
filen  = pthfilen.split('/')[-1]
plotdir = 'frames'
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

h5f = h5py.File(pthfilen,'r')
if len(sys.argv) < 3:
    print("Available datafields: ", list(h5f['field_types'].keys()))
    exit(1)
if sys.argv[2] == "_all_":
    varlist = h5f['field_types'].keys()
else:
    varlist  = sys.argv[2:]

print(varlist)

print("Reading file: %s" % pthfilen)
for var in varlist:
    #output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
    fnl = filen.split('/')[-1]
    output = plotdir+'/'+'_'.join(fnl.split('_')[:-1])+'_'+var+'_'+fnl.split('_')[-1].replace('.h5',".png")
    pc.plotcompose(pthfilen, var, output)

h5f.close()
