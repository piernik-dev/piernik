#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import h5py
from colored_io import die, prtinfo, prtwarn

# This script reads problem.par field from provided
# h5/res file, extracts it and saves it in CWD,
# so it can be used in new run.

if (len(sys.argv) > 1):
    hdf5_filename = sys.argv[1]
    prtinfo("File to process is %s" % hdf5_filename)
else:
    die("No file provided, exiting.")

try:
    fh5 = h5py.File(hdf5_filename, 'r')
    parfile_out = open("problem.par.%s" % hdf5_filename.strip(".res").strip(".h5"), 'w')
except:
    die("Problem opening %s file. " % hdf5_filename)

parameter_object = fh5["problem.par"]
for line in parameter_object:
    # print(line.decode("utf-8"))
    parfile_out.write(line.decode("utf-8") + "\n")

prtinfo("Successfully read parameters from %s file, list of parameters saved to %s" % (hdf5_filename, parfile_out.name))
parfile_out.close()
fh5.close()
